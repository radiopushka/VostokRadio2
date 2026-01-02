#include "FFT.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct FFT_rsmp *FFT_resample_init(int bins,int OF_bins,int ring_buffer_delay, float fs, float fend,float srate){
    struct FFT_rsmp *rsmp = malloc(sizeof(struct FFT_rsmp));
    float freq_pbin = (fend-fs)/((float)bins);
    rsmp->alpha = freq_pbin/srate;
    rsmp->nalpha = 1-rsmp->alpha;
    rsmp->bins = bins;

    rsmp->OF_scale = 1.0f/bins;



    int fft_buff_len = srate+240;//prevent clicking so that all synthesized signals are being sampled

    rsmp->hpf_alpha = (10.0f)/srate;//10hz high pass
    rsmp->hpf_nalpha = 1-rsmp->hpf_alpha;

    //how many bins starting from the top of the band are resampled orthogonally to the MPX signals
    //asumes the 19khz pilot tone is sampled as a sine function
    rsmp->OFbins = OF_bins<<1;


    rsmp->ring_buffer_size = ring_buffer_delay;
    rsmp->rastoyanee = malloc(sizeof(int)*bins);
    memset(rsmp->rastoyanee, 0, sizeof(int)*bins);
    rsmp->count_buff = malloc(sizeof(int)*ring_buffer_delay);
    memset(rsmp->count_buff, 0, sizeof(int)*ring_buffer_delay);
    rsmp->count_buff_i = rsmp->count_buff;
    rsmp->count_buff_e = rsmp->count_buff + ring_buffer_delay;


    rsmp->pre_high_pass = 0;
    rsmp->lpf_r = malloc(bins*sizeof(float));
    rsmp->lpf_i = malloc(bins*sizeof(float));
    memset(rsmp->lpf_r, 0, bins*sizeof(float));
    memset(rsmp->lpf_i, 0, bins*sizeof(float));
    rsmp->lpf_r_v = malloc(ring_buffer_delay*sizeof(float*));
    rsmp->lpf_i_v = malloc(ring_buffer_delay*sizeof(float*));
    for(int i = 0;i<ring_buffer_delay;i++){
        rsmp->lpf_r_v[i] = malloc(sizeof(float)*bins);
        rsmp->lpf_i_v[i] = malloc(sizeof(float)*bins);
        memset(rsmp->lpf_r_v[i], 0, bins*sizeof(float));
        memset(rsmp->lpf_i_v[i], 0, bins*sizeof(float));
    }
    rsmp->lpf_i_l = rsmp->lpf_i_v + ring_buffer_delay;
    rsmp->lpf_r_l = rsmp->lpf_r_v + ring_buffer_delay;
    rsmp->lpf_i_i = rsmp->lpf_i_v;
    rsmp->lpf_r_i = rsmp->lpf_r_v;

    rsmp->amplitude = malloc(bins*sizeof(float));
    memset(rsmp->amplitude, 0, bins*sizeof(float));

    rsmp->amplitude_ring = malloc(ring_buffer_delay*sizeof(float*));
    for(int i = 0;i<ring_buffer_delay;i++){
        rsmp->amplitude_ring[i] = malloc(sizeof(float)*bins);
        memset(rsmp->amplitude_ring[i], 0, bins*sizeof(float));
    }
    rsmp->amplitude_ring_i = rsmp->amplitude_ring;
    rsmp->amplitude_ring_end = rsmp->amplitude_ring + ring_buffer_delay;




    rsmp->length = fft_buff_len;
    rsmp->counter = 0;

    int size_dsp = (rsmp->length*bins)*2;
    rsmp->dsp_array = malloc(size_dsp*sizeof(float));
    /*for(float* i = rsmp->dsp_array;i<rsmp->dsp_array+size_dsp;i++){
        *i = 0;

    }*/

    float fcnt = fs;
    for(int i = 0;i<bins;i++){
        float shifter = (fcnt/(srate * 4))*(2*M_PI);

        long double counter = 0;
        for(int i2 = 0;i2<rsmp->length;i2++){

            //resample by 4
            long double ac = 0;
            long double as = 0;
            for(int i3 = 0;i3<4;i3++){
                ac = ac + cosl(counter);
                as = as + sinl(counter);
                counter += shifter;
                if(counter > M_PI*2){
                    counter -= M_PI*2;
                }
            }
            ac = ac/4.0;
            as = as/4.0;

            int location = i2*(bins*2) + i*2;
            float* cval = rsmp->dsp_array + location;
            float* sval = rsmp->dsp_array + location + 1;
            *cval = ac;
            *sval = as;
        }


        fcnt = fcnt + freq_pbin;
    }



    return rsmp;
}

int get_resamp_size(struct FFT_rsmp *rsmp){
    return rsmp->bins;
}
void calc_max_buffer_amp(struct FFT_rsmp *rsmp){
    //пересчитываем максимальную амплитуду из буфера
    for(int* restrict  i = rsmp->rastoyanee;i<rsmp->rastoyanee+rsmp->bins;i++){
        *i = rsmp->ring_buffer_size;
    }
    float* a_end = rsmp->amplitude + rsmp->bins;
    float** restrict i_camp = rsmp->amplitude_ring_i+1;
    for(int i = 1;i<rsmp->ring_buffer_size;i++,i_camp++){
        if(i_camp >= rsmp->amplitude_ring_end)
            i_camp = rsmp->amplitude_ring;
        float* restrict amp_list = *i_camp;
        int* restrict rastoyanee = rsmp->rastoyanee;
        for(float* restrict i2 = rsmp->amplitude;i2<a_end;i2++,amp_list++,rastoyanee++){
            float lamp = *amp_list;
            if(lamp > *i2){
                *rastoyanee = (rsmp->ring_buffer_size) - i;
                *i2 = lamp;
            }
        }
    }
}
float* resamp_pre_process(struct FFT_rsmp *rsmp, float in,float* restrict eq){
    float* restrict lpr = rsmp->lpf_r;
    float* restrict lpi = rsmp->lpf_i;
    float* restrict amp = rsmp->amplitude;

    float* restrict darray = rsmp->dsp_array;

    int cnt = rsmp->counter;
    cnt++;
    if(cnt >= rsmp->length){
     cnt = 0;
    }
    int index = cnt*((rsmp->bins)<<1);
    rsmp->counter = cnt;


    rsmp->pre_high_pass = rsmp->pre_high_pass*rsmp->hpf_nalpha + in*rsmp->hpf_alpha;

    float mod = in - rsmp->pre_high_pass;

    float* dend = darray + index + (rsmp->bins<<1);
    for(float* restrict ic = darray + index;ic<dend;ic = ic+2,lpi++,lpr++,amp++,eq++){
        float r = mod*(*ic);
        float i = mod*(*(ic+1));
        r = r*(*eq);
        i = i*(*eq);

        *lpi = (*lpi)*rsmp->nalpha + i*(rsmp->alpha);
        *lpr = (*lpr)*rsmp->nalpha + r*(rsmp->alpha);
        //speed optimization
        *amp = sqrtf((float)((*lpi)*(*lpi) + (*lpr)*(*lpr)));

    }

    //sdvig_vverh(rsmp);

    rsmp->lpf_r_i = rsmp->lpf_r_i + 1;
    rsmp->lpf_i_i = rsmp->lpf_i_i + 1;
    if(rsmp->lpf_r_i >= rsmp->lpf_r_l)
        rsmp->lpf_r_i = rsmp->lpf_r_v;
    if(rsmp->lpf_i_i >= rsmp->lpf_i_l)
        rsmp->lpf_i_i = rsmp->lpf_i_v;
    memcpy(*(rsmp->lpf_r_i),rsmp->lpf_r,rsmp->bins*sizeof(float));
    memcpy(*(rsmp->lpf_i_i),rsmp->lpf_i,rsmp->bins*sizeof(float));


    rsmp->amplitude_ring_i = rsmp->amplitude_ring_i + 1;
    if(rsmp->amplitude_ring_i >= rsmp->amplitude_ring_end)
        rsmp->amplitude_ring_i = rsmp->amplitude_ring;
    memcpy(*(rsmp->amplitude_ring_i),rsmp->amplitude,rsmp->bins*sizeof(float));
    //now find the maximum amplitude for each frequency and calculate the number of samples before it reaches get_signal

    rsmp->count_buff_i = rsmp->count_buff_i + 1;
    if(rsmp->count_buff_i >= rsmp->count_buff_e)
        rsmp->count_buff_i = rsmp->count_buff;
    *(rsmp->count_buff_i)=index;

    calc_max_buffer_amp(rsmp);
    return rsmp->amplitude;
}

float resamp_get_signal(struct FFT_rsmp *rsmp, float* eq){
    float** restrict lp_r = rsmp->lpf_r_i+1;
    if(lp_r >= rsmp->lpf_r_l)
        lp_r = rsmp->lpf_r_v;
    float* restrict lpr = *lp_r;

    float** restrict lp_i = rsmp->lpf_i_i+1;
    if(lp_i >= rsmp->lpf_i_l)
        lp_i = rsmp->lpf_i_v;
    float* restrict lpi = *lp_i;


    float* restrict darray = rsmp->dsp_array;

    int* count = rsmp->count_buff_i+1;
    if(count >= rsmp->count_buff_e)
        count = rsmp->count_buff;


    int cnt = *count;

    float* dend = darray + cnt + (rsmp->bins<<1);
    //Orthogonal highs
    float* last_bin = dend-rsmp->OFbins;

    float OF_scale = 1.0f - rsmp->OF_scale;
    float output = 0.0f;
    for(float* restrict ic = darray + cnt;ic<dend;ic = ic+2,lpi++,lpr++,eq++){


        float sval;
        if(last_bin <= ic){
            // the 19 khz pilot is generated by a sine function
            // only resample the cosine part of the fourie image
            // reduce the imaginary part
            // should help but far from ideal
            sval = (*lpr)*(*ic) + (*lpi)*(*(ic+1))*OF_scale;
            if(OF_scale > 0.0f)
                OF_scale = OF_scale - rsmp->OF_scale;
            else
                OF_scale = 0.0f;
        }else{
            sval = (*lpi)*(*(ic+1)) + (*lpr)*(*ic);
        }
        output += sval*(*eq);
    }
    return output;
}

void free_resamp(struct FFT_rsmp *rsmp){
    free(rsmp->lpf_r);
    free(rsmp->lpf_i);
    free(rsmp->amplitude);
    free(rsmp->dsp_array);
    for(int i = 0;i<rsmp->ring_buffer_size;i++){
        free(rsmp->lpf_r_v[i]);
        free(rsmp->lpf_i_v[i]);
        free(rsmp->amplitude_ring[i]);
    }
    free(rsmp->count_buff);
    free(rsmp->rastoyanee);
    free(rsmp->lpf_r_v);
    free(rsmp->lpf_i_v);
    free(rsmp->amplitude_ring);
    free(rsmp);

}
