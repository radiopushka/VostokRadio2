#include "alsa_pipe/alsa_pipe.h"
#include "./FFT/FFT.h"
#include "./EQ_max_algorithm/maxim.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

//the Raspberry Pi will have no cli control, purely config file based
//CPU is limited
#include "./config_file/config.h"


char recording[32];
char playback[32];

const float int_value = 2147483647.0;
//rates fixed to 48khz and 192 khz
//
//user settings:
float pilot_amp = 0.15;
int OF_bins = 1;
//post AGC limiter amp and image spectral limiter
float pre_amp = 0.7;
float post_amp =  1;
float limit_value = 2.3e9;
float* l_release;
float harmonic_diff = 0.6;
int lookahead = 3;
//the raspberry pi zero can do 9 bins max
int bins = 9;//valid values: 5,9,15,30,45
//pre limiter equalization
float* pre_eq;
//gain controller
float attack = 0.03;
float release = 0.003;
float target = 6e9;//3e9 for pi zero and 6e9 for normal setups
float noise_th = 2e6;


//mono to st_ratio
float stereo_ratio = 0.3;


//high pass filter
float alpha = (400.0)/48000.0;
float nalpha;//10hz
float hpv_l = 0;
float hpv_r = 0;

//bass boost
float bhpv_l = 0;
float bhpv_r = 0;
float bass_boost = 0.4;
float nbass_boost;
//so that we can call this later from a gui
void pre_set_settings(){
    nalpha = 1-alpha;
    nbass_boost = 1-bass_boost;
}
void setup_globals_from_config(char* file){
    CfgRaster cfg = read_cfg_file(file);
    if(!cfg){
        printf("config file not found\n");
        return;
    }

    //audio playback recording
    char* recording_f = get_value_by(cfg,"audio","iface");
    char* playback_f = get_value_by(cfg,"audio","oface");
    if(recording_f)
        sprintf(recording,"%s",recording_f);
    if(playback_f)
        sprintf(playback,"%s",playback_f);
    printf("read config\n");
    printf("playback: %s, rec: %s\n",playback,recording);

    //MPX settings
    char* pilot_amp_f = get_value_by(cfg,"MPX","pilot_amp");
    if(pilot_amp_f){
        pilot_amp = atof(pilot_amp_f);
        printf("pilot amp: %f\n",pilot_amp);
    }

    char* stereo_ratio_f = get_value_by(cfg,"MPX","stereo_ratio");
    if(stereo_ratio_f){
        stereo_ratio = atof(stereo_ratio_f);
        printf("stereo ratio: %f\n",stereo_ratio);
    }
    char* OF_bins_f = get_value_by(cfg,"MPX","OF_bins");
    if(OF_bins_f){
        OF_bins = atof(OF_bins_f);
        printf("Orthogonal(to 19 khz pilot) Frequency bins: %d\n",OF_bins);
    }


    //limiter
    char* pre_amp_f = get_value_by(cfg,"limiter","pre_amp");
    if(pre_amp_f){
        pre_amp = atof(pre_amp_f);
        printf("pre amp: %f\n",pre_amp);
    }
    char* lookahead_f = get_value_by(cfg,"limiter","lookahead");
    if(lookahead_f){
        lookahead = atoi(lookahead_f);
        printf("lookahead: %d\n",lookahead);
    }
    char* post_amp_f = get_value_by(cfg,"limiter","post_amp");
    if(post_amp_f){
        post_amp = atof(post_amp_f);
        printf("post amp: %f\n",post_amp);
    }
    char* limit_value_f = get_value_by(cfg,"limiter","limit");
    if(limit_value_f){
        limit_value = atof(limit_value_f);
        printf("limit: %f\n",limit_value);
    }
    char* harmonic_diff_f = get_value_by(cfg,"limiter","save_harmonics");
    if(harmonic_diff_f){
        harmonic_diff = atof(harmonic_diff_f);
        printf("harmonic diff: %f\n",harmonic_diff);
    }

    char* bins_f = get_value_by(cfg,"limiter","bins");
    if(bins_f){
        int bin_r = atoi(bins_f);
        if(bin_r == 5 || bin_r == 9 || bin_r == 15 || bin_r == 30 || bin_r == 45){
            bins = bin_r;
        }
        printf("bins: %d\n",bins);
    }

    //gain controller
    char* attack_f = get_value_by(cfg,"agc","attack");
    if(attack_f){
        attack = atof(attack_f);
        printf("attack: %f\n",attack);
    }
    char* release_f = get_value_by(cfg,"agc","release");
    if(release_f){
        release = atof(release_f);
        printf("release: %f\n",release);
    }
    char* target_f = get_value_by(cfg,"agc","target");
    if(target_f){
        target = atof(target_f);
        printf("target: %f\n",target);
    }
    char* noise_th_f = get_value_by(cfg,"agc","noise_th");
    if(noise_th_f){
        noise_th = atof(noise_th_f);
        printf("noise th: %f\n",noise_th);
    }

    //high pass filter
    char* alpha_f = get_value_by(cfg,"highpass","alpha");
    if(alpha_f){
        alpha = atof(alpha_f);
        printf("alpha: %f\n",alpha);
    }
    //bass boost
    char* bass_boost_f = get_value_by(cfg,"bassboost","boost");
    if(bass_boost_f){
        bass_boost = atof(bass_boost_f);
        printf("bass boost: %f\n",bass_boost);
    }

    free_cfg_mem(cfg);

}
//setup the eq and limiter release
void setup_from_config(float* eq,float* release,char* file){

    CfgRaster cfg = read_cfg_file(file);
    if(!cfg)
        return;
    char val[13];
    for(int i = 0;i<bins;i++){
        sprintf(val,"%d",i);

        char* e_val = get_value_by(cfg,"eq",val);
        char* r_val = get_value_by(cfg,"release",val);
        if(e_val){
            eq[i] = atof(e_val);
        }
        if(r_val){
            release[i] = atof(r_val);
        }

    }

    free_cfg_mem(cfg);

}
int main(int argn,char* argv[]){

    if(argn >= 2){
            setup_globals_from_config(argv[1]);
    }else{
        sprintf(recording,"default");
        sprintf(playback,"default");
    }
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
    float* midbuff_m = malloc(sizeof(float)*half_b);
    float* midbuff_s = malloc(sizeof(float)*half_b);
    memset(midbuff_m,0,sizeof(float)*half_b);
    memset(midbuff_s,0,sizeof(float)*half_b);
    int* output = malloc(sizeof(int)*rate1);


    //equalizer, programatically controlled
    float eq[bins];
    float eq_s[bins];
    float eq_helper[bins];
    float eq_helper_s[bins];
    pre_eq = malloc(sizeof(float)*bins);
    l_release = malloc(sizeof(float)*bins);
    for(int i = 0;i<bins;i++){
        eq[i] = 1;
        eq_s[i] = 1;
        pre_eq[i] = 1;
        l_release[i] = 0.5;
    }
    //user setting eq

    //default setting eq 45 bins
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

    //default setting for eq 9 bins
    if(bins == 9){
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
    }
    //user setting eq and release
    if(argn >= 2){
        setup_from_config(eq,l_release,argv[1]);
    }

    //FFT resampling mono
    struct FFT_rsmp *rsmp = FFT_resample_init(bins,OF_bins,lookahead, 1000, 16000, rate1);
    struct FFT_rsmp *rsmp_st = FFT_resample_init(bins,OF_bins,lookahead, 1000, 16000, rate1);
    //gain controller
    struct Gain_Control *gc = gain_control_init(attack,release,target,noise_th);



    //mpx syntheizer

    long double shifter_19 = ((19000.0) / (rate2*4.0))*(2*M_PI);
    long double shifter_38 = shifter_19*2;

    float* synth_19 = malloc(sizeof(float)*mpx_b);
    float* synth_38 = malloc(sizeof(float)*mpx_b);

    float* synth_38_i = synth_38;
    float* synth_19_i = synth_19;
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
        float* i_mb = midbuff_m;
        float* i_sb = midbuff_s;
        set_gain_control(gc,target,attack,release,noise_th);
        for(int* sp = recbuff;sp<recbuff_end;sp = sp+2){

            float l = ((float)(*sp))*pre_amp;
            float r = ((float)(*(sp+1)))*pre_amp;
            hpv_r = hpv_r*nalpha+r*alpha;
            hpv_l = hpv_l*nalpha+l*alpha;
            r = r - hpv_r;
            l = l - hpv_l;
            gain_control(gc,&l,&r);
            bhpv_r = bhpv_r*nalpha+r*alpha;
            bhpv_l = bhpv_l*nalpha+l*alpha;
            l = l*nbass_boost+bhpv_l*bass_boost;
            r = r*nbass_boost+bhpv_r*bass_boost;
            float sum = l+r;
            float diff = l-r;
            *i_mb = sum;
            *i_sb = diff;
            i_mb++;
            i_sb++;

        }
        //calculate the i and r FFT coefficients for each
        i_mb = midbuff_m;
        i_sb = midbuff_s;
        float suma;

        for(int i = 0;i<half_b;i++){

            float* amplitude = resamp_pre_process(rsmp, *i_mb,pre_eq);
            suma = find_amp(amplitude,bins);

            adjust_eq(eq,amplitude,rsmp->rastoyanee,bins,limit_value,suma,eq_helper,l_release,harmonic_diff);
            float val = resamp_get_signal(rsmp,eq);


            amplitude = resamp_pre_process(rsmp_st, *i_sb,pre_eq);
            suma = find_amp(amplitude,bins);

            adjust_eq(eq_s,amplitude,rsmp_st->rastoyanee,bins,limit_value,suma,eq_helper_s,l_release,harmonic_diff);
            float stval = resamp_get_signal(rsmp_st,eq_s);

            *i_mb = val*post_amp;
            *i_sb = stval*post_amp;
            i_mb++;
            i_sb++;

        }
        //sample up to 192khz and add the MPX signals

        int* oi = output;
        i_mb = midbuff_m;
        i_sb = midbuff_s;

        float limit_audio = (1-pilot_amp);
        float pilot_v = pilot_amp*int_value;
        for(int i = 0;i<half_b;i++){
            float mono = *i_mb;
            float stereo = *i_sb;
            i_mb++;i_sb++;

            float stereo_val = stereo*stereo_ratio;
            float headroom = int_value - fabs(stereo_val);
            float ratio_mono = headroom/int_value;
            mono = mono*ratio_mono;
            stereo = stereo_val;



            //put 8 samples for 192 khz
            for(int i2 = 0;i2<4;i2++){
                float w38 = synth_38[mpx_count];
                float w19 = synth_19[mpx_count];

                float sample = w19*pilot_v+(w38*stereo+mono)*limit_audio;
                //float sample = stereo+mono;
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
