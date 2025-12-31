#include "FFT.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct FFT_rsmp *FFT_resample_init(int bins,int ring_buffer_delay, float fs, float fend,float srate){
    struct FFT_rsmp *rsmp = malloc(sizeof(struct FFT_rsmp));
    float freq_pbin = (fend-fs)/((float)bins);
    rsmp->alpha = freq_pbin/srate;
    rsmp->nalpha = 1-rsmp->alpha;
    rsmp->bins = bins;


    int fft_buff_len = srate+240;//prevent clicking so that all synthesized signals are being sampled

    rsmp->hpf_alpha = (10.0f)/srate;//10hz high pass
    rsmp->hpf_nalpha = 1-rsmp->hpf_alpha;


    rsmp->ring_buffer_size = ring_buffer_delay;
    rsmp->rastoyanee = malloc(sizeof(int)*bins);
    memset(rsmp->rastoyanee, 0, sizeof(int)*bins);
    rsmp->count_buff = malloc(sizeof(int)*ring_buffer_delay);
    memset(rsmp->count_buff, 0, sizeof(int)*ring_buffer_delay);


    rsmp->pre_high_pass = 0;
    rsmp->lpf_r = malloc(bins*sizeof(double));
    rsmp->lpf_i = malloc(bins*sizeof(double));
    rsmp->lpf_r_v = malloc(ring_buffer_delay*sizeof(double*));
    rsmp->lpf_i_v = malloc(ring_buffer_delay*sizeof(double*));
    for(int i = 0;i<ring_buffer_delay;i++){
        rsmp->lpf_r_v[i] = malloc(sizeof(double)*bins);
        rsmp->lpf_i_v[i] = malloc(sizeof(double)*bins);
        memset(rsmp->lpf_r_v[i], 0, bins*sizeof(double));
        memset(rsmp->lpf_i_v[i], 0, bins*sizeof(double));
    }

    rsmp->amplitude = malloc(bins*sizeof(double));
    memset(rsmp->amplitude, 0, bins*sizeof(double));

    rsmp->amplitude_ring = malloc(ring_buffer_delay*sizeof(double*));
    for(int i = 0;i<ring_buffer_delay;i++){
        rsmp->amplitude_ring[i] = malloc(sizeof(double)*bins);
        memset(rsmp->amplitude_ring[i], 0, bins*sizeof(double));
    }



    rsmp->length = fft_buff_len;
    rsmp->counter = 0;

    rsmp->cos_array = malloc(bins*sizeof(double*));
    rsmp->sin_array = malloc(bins*sizeof(double*));

    float fcnt = fs;
    for(int i = 0;i<bins;i++){
        double shifter = (fcnt/(srate * 4))*(2*M_PI);
        rsmp->cos_array[i] = malloc(sizeof(double)*rsmp->length);
        rsmp->sin_array[i] = malloc(sizeof(double)*rsmp->length);

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
            rsmp->cos_array[i][i2] = ac;
            rsmp->sin_array[i][i2] = as;
        }


        fcnt = fcnt + freq_pbin;
    }



    return rsmp;
}

int get_resamp_size(struct FFT_rsmp *rsmp){
    return rsmp->bins;
}
void sdvig_vverh(struct FFT_rsmp *rsmp){
    //сдвигаем массивы с FFT коэфф в конец на один шаг

    for(int s = rsmp->ring_buffer_size-1;s>0;s--){
        memmove(rsmp->lpf_r_v[s],rsmp->lpf_r_v[s-1],rsmp->bins*sizeof(double));
        memmove(rsmp->lpf_i_v[s],rsmp->lpf_i_v[s-1],rsmp->bins*sizeof(double));
        memmove(rsmp->amplitude_ring[s],rsmp->amplitude_ring[s-1],rsmp->bins*sizeof(double));
        rsmp->count_buff[s] = rsmp->count_buff[s-1];
    }
}
void calc_max_buffer_amp(struct FFT_rsmp *rsmp){
    //пересчитываем максимальную амплитуду из буфера
    for(int i = 0;i<rsmp->bins;i++){
        rsmp->rastoyanee[i] = rsmp->ring_buffer_size;
    }
    for(int i = 1;i<rsmp->ring_buffer_size;i++){
        double* amp_list = rsmp->amplitude_ring[i];
        for(int i2 = 0;i2<rsmp->bins;i2++){
            double lamp = amp_list[i2];
            if(lamp > rsmp->amplitude[i2]){
                rsmp->rastoyanee[i2] = (rsmp->ring_buffer_size) - i;
                rsmp->amplitude[i2] = lamp;
            }
        }
    }
}
double* resamp_pre_process(struct FFT_rsmp *rsmp, double in,double* restrict eq){
    double* restrict lpr = rsmp->lpf_r;
    double* restrict lpi = rsmp->lpf_i;
    double* restrict amp = rsmp->amplitude;

    double** restrict sarray = rsmp->sin_array;
    double** restrict carray = rsmp->cos_array;

    int cnt = rsmp->counter;
    cnt++;
    if(cnt >= rsmp->length){
     cnt = 0;
    }
    rsmp->counter = cnt;

    rsmp->pre_high_pass = rsmp->pre_high_pass*rsmp->hpf_nalpha + in*rsmp->hpf_alpha;

    double mod = in - rsmp->pre_high_pass;

    for(int ic = 0;ic<rsmp->bins;ic++){
        double r = mod*carray[ic][cnt];
        double i = mod*sarray[ic][cnt];
        r = r*(*eq);
        i = i*(*eq);

        *lpi = (*lpi)*rsmp->nalpha + i*(rsmp->alpha);
        *lpr = (*lpr)*rsmp->nalpha + r*(rsmp->alpha);
        *amp = sqrt((*lpi)*(*lpi) + (*lpr)*(*lpr));
        lpi++; lpr++; amp++;
        eq++;

    }

    sdvig_vverh(rsmp);
    memmove(rsmp->lpf_r_v[0],rsmp->lpf_r,rsmp->bins*sizeof(double));
    memmove(rsmp->lpf_i_v[0],rsmp->lpf_i,rsmp->bins*sizeof(double));
    memmove(rsmp->amplitude_ring[0],rsmp->amplitude,rsmp->bins*sizeof(double));
    //now find the maximum amplitude for each frequency and calculate the number of samples before it reaches get_signal
    rsmp->count_buff[0] = cnt;
    calc_max_buffer_amp(rsmp);
    return rsmp->amplitude;
}

double resamp_get_signal(struct FFT_rsmp *rsmp, double* eq){
    double* restrict lpr = rsmp->lpf_r_v[rsmp->ring_buffer_size-1];
    double* restrict lpi = rsmp->lpf_i_v[rsmp->ring_buffer_size-1];

    double** restrict sarray = rsmp->sin_array;
    double** restrict carray = rsmp->cos_array;

    int cnt = rsmp->count_buff[rsmp->ring_buffer_size-1];

    double output = 0;
    for(int ic = 0;ic<rsmp->bins;ic++){

        double sval = (*lpi)*sarray[ic][cnt] + (*lpr)*carray[ic][cnt];
        lpi++;
        lpr++;
        output += sval*(*eq);
        eq++;
    }
    return output;
}

void free_resamp(struct FFT_rsmp *rsmp){
    free(rsmp->lpf_r);
    free(rsmp->lpf_i);
    free(rsmp->amplitude);
    for(int i = 0;i<rsmp->bins;i++){
        free(rsmp->cos_array[i]);
        free(rsmp->sin_array[i]);
    }
    for(int i = 0;i<rsmp->ring_buffer_size;i++){
        free(rsmp->lpf_r_v[i]);
        free(rsmp->lpf_i_v[i]);
        free(rsmp->amplitude_ring[i]);
    }
    free(rsmp->count_buff);
    free(rsmp->rastoyanee);
    free(rsmp->lpf_r_v);
    free(rsmp->lpf_r_v);
    free(rsmp->lpf_i_v);
    free(rsmp->amplitude_ring);
    free(rsmp->cos_array);
    free(rsmp->sin_array);

}
