#include "alsa_pipe/alsa_pipe.h"
#include "./FFT/FFT.h"
#include "./EQ_max_algorithm/maxim.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>


char* recording = "default";
char* playback ="default";

const double int_value = 2147483647.0;
//rates fixed to 48khz and 192 khz
//
//user settings:
double pilot_amp = 0.15;
//post AGC limiter amp and image spectral limiter
double pre_amp = 0.7;
double post_amp =  1;
double limit_value = 2.3e9;
double* l_release;
double harmonic_diff = 0.6;

//the raspberry pi zero can do 9 bins max
int bins = 9;//valid values: 5,9,15,30,45
//pre limiter equalization
double* pre_eq;
//gain controller
double attack = 0.01;
double release = 0.001;
double target = 6e9;//3e9 for pi zero and 6e9 for normal setups
double noise_th = 2e6;


//mono to st_ratio
double stereo_ratio = 0.3;


//high pass filter
double alpha = (400.0)/48000.0;
double nalpha;//10hz
double hpv_l = 0;
double hpv_r = 0;

//bass boost
double bhpv_l = 0;
double bhpv_r = 0;
double bass_boost = 0.4;
double nbass_boost;
//so that we can call this later from a gui
void pre_set_settings(){
    nalpha = 1-alpha;
    nbass_boost = 1-bass_boost;
}
int main(){
    pre_set_settings();

    int mpx_b = 189999;

    int ch1 = 2;
    int ch2 = 2;
    int rate1 = 48000;// mandatory
    int rate2 = 192000;// mandatory
    int input_buffer_prop = rate2/rate1;
    int i_buffer_size = rate1/input_buffer_prop;
    int half_b = i_buffer_size/2;

    int* recbuff = malloc(sizeof(int)*i_buffer_size);
    int* recbuff_end = recbuff+i_buffer_size;
    memset(recbuff,0, sizeof(int)*i_buffer_size);
    double* midbuff_m = malloc(sizeof(double)*half_b);
    double* midbuff_s = malloc(sizeof(double)*half_b);
    memset(midbuff_m,0,sizeof(double)*half_b);
    memset(midbuff_s,0,sizeof(double)*half_b);
    int* output = malloc(sizeof(int)*rate1);


    //equalizer, programatically controlled
    double eq[bins];
    double eq_s[bins];
    double eq_helper[bins];
    double eq_helper_s[bins];
    pre_eq = malloc(sizeof(double)*bins);
    l_release = malloc(sizeof(double)*bins);
    for(int i = 0;i<bins;i++){
        eq[i] = 1;
        eq_s[i] = 1;
        pre_eq[i] = 1;
        l_release[i] = 0.5;
    }
    //user setting eq for 45 bins
    /*
    pre_eq[0]=0.0001;
    pre_eq[1]=0.001;
    pre_eq[2]=0.001;
    pre_eq[3]=0.001;
    pre_eq[4]=0.001;
    pre_eq[5]=0.1;
    pre_eq[6]=0.6;
    pre_eq[7]=0.8;
    pre_eq[8]=0.9;

    pre_eq[40]=1.5;
    pre_eq[41]=1.5;
    pre_eq[43]=1.5;
    pre_eq[44]=1.5;
    */

    //user setting for eq 9 bins
    pre_eq[0]=0.005;
    pre_eq[1]=0.2;
    pre_eq[2]=0.7;
    pre_eq[3]=0.75;
    pre_eq[4]=0.78;
    pre_eq[5]=0.8;
    pre_eq[6]=1;
    pre_eq[7]=1.1;
    pre_eq[8]=1.1;
    //release settings:
    l_release[0]=0.0000001;
    l_release[1]=0.0001;

    //FFT resampling mono
    struct FFT_rsmp *rsmp = FFT_resample_init(bins,3, 1000, 16000, rate1);
    struct FFT_rsmp *rsmp_st = FFT_resample_init(bins,3, 1000, 16000, rate1);
    //gain controller
    struct Gain_Control *gc = gain_control_init(attack,release,target,noise_th);



    //mpx syntheizer

    long double shifter_19 = ((19000.0) / (rate2*4.0))*(2*M_PI);
    long double shifter_38 = shifter_19*2;

    double* synth_19 = malloc(sizeof(double)*mpx_b);
    double* synth_38 = malloc(sizeof(double)*mpx_b);

    double* synth_38_i = synth_38;
    double* synth_19_i = synth_19;
    long double count = 0;
    for(int i = 0;i<mpx_b;i++){
        long double a19 = 0;
        long double a38 = 0;
        for(int i = 0;i<4;i++){
            a19 = a19 + sinl(shifter_19*count);
            a38 = a38 + sinl(shifter_38*count);
            count = count + 1.0;
        }
        *synth_19_i = a19/4.0;
        *synth_38_i = a38/4.0;
        synth_19_i++;
        synth_38_i++;

    }


    if(setup_alsa_pipe(recording, playback, &ch1, &ch2,&rate1, &rate2, rate1)== -1){
        printf("alsa error\n");
        goto exit;
    }else{
        printf("started main thread\n");
    }



    int mpx_count = 0;
    while(1){

        if(get_audio(recbuff,i_buffer_size)== -1)
            break;

        //break the stereo signal into the L+R and L-R buffers
        double* i_mb = midbuff_m;
        double* i_sb = midbuff_s;
        set_gain_control(gc,target,attack,release,noise_th);
        for(int* sp = recbuff;sp<recbuff_end;sp = sp+2){

            double l = ((double)(*sp))*pre_amp;
            double r = ((double)(*(sp+1)))*pre_amp;
            hpv_r = hpv_r*nalpha+r*alpha;
            hpv_l = hpv_l*nalpha+l*alpha;
            r = r - hpv_r;
            l = l - hpv_l;
            gain_control(gc,&l,&r);
            bhpv_r = bhpv_r*nalpha+r*alpha;
            bhpv_l = bhpv_l*nalpha+l*alpha;
            l = l*nbass_boost+bhpv_l*bass_boost;
            r = r*nbass_boost+bhpv_r*bass_boost;
            double sum = l+r;
            double diff = l-r;
            *i_mb = sum;
            *i_sb = diff;
            i_mb++;
            i_sb++;

        }
        //calculate the i and r FFT coefficients for each
        i_mb = midbuff_m;
        i_sb = midbuff_s;
        double suma;

        for(int i = 0;i<half_b;i++){

            double* amplitude = resamp_pre_process(rsmp, *i_mb,pre_eq);
            suma = find_amp(amplitude,bins);

            adjust_eq(eq,amplitude,rsmp->rastoyanee,bins,limit_value,suma,eq_helper,l_release,harmonic_diff);
            double val = resamp_get_signal(rsmp,eq);


            amplitude = resamp_pre_process(rsmp_st, *i_sb,pre_eq);
            suma = find_amp(amplitude,bins);

            adjust_eq(eq_s,amplitude,rsmp_st->rastoyanee,bins,limit_value,suma,eq_helper_s,l_release,harmonic_diff);
            double stval = resamp_get_signal(rsmp_st,eq_s);

            *i_mb = val*post_amp;
            *i_sb = stval*post_amp;
            i_mb++;
            i_sb++;

        }
        //sample up to 192khz and add the MPX signals

        int* oi = output;
        i_mb = midbuff_m;
        i_sb = midbuff_s;

        double limit_audio = (1-pilot_amp);
        double pilot_v = pilot_amp*int_value;
        for(int i = 0;i<half_b;i++){
            double mono = *i_mb;
            double stereo = *i_sb;
            i_mb++;i_sb++;

            double stereo_val = stereo*stereo_ratio;
            double headroom = int_value - fabs(stereo_val);
            double ratio_mono = headroom/int_value;
            mono = mono*ratio_mono;
            stereo = stereo_val;



            //put 8 samples for 192 khz
            for(int i2 = 0;i2<4;i2++){
                double w38 = synth_38[mpx_count];
                double w19 = synth_19[mpx_count];

                double sample = w19*pilot_v+(w38*stereo+mono)*limit_audio;
                //double sample = stereo+mono;
                *oi = sample;
                oi++;
                *oi = sample;
                oi++;
                mpx_count++;
                if(mpx_count >= mpx_b)
                    mpx_count = 0;
            }
        }
        queue_audio(output);
    }

exit:
    free_resamp(rsmp);
    free_resamp(rsmp_st);

    free_gain_control(gc);
    free(pre_eq);
    free(recbuff);
    free(midbuff_m);
    free(midbuff_s);
    free(output);
    free(synth_19);
    free(synth_38);
    alsa_pipe_exit();
}
