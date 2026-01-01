#include "maxim.h"
#include <stdlib.h>
#include <math.h>
#include <string.h>


float find_amp(float *fft_out, int bins){
    float suma = 0;
    for(float* restrict ittr = fft_out; ittr < fft_out + bins; ittr++){
        suma = suma + (*ittr);
    }
    return suma;
}
//eq ratio compression:
/*
 * compression = 0.1 //the lower this number the less the amplitude difference
 * elements = bins //number of elements in eq
 * //for each frequency:
 * f1_ratio = f1_ratio_old * compression + (1 - compression) / elements
 *
 */
void adjust_eq(float *eq,float* fft_out,int* rastoyane, int bins, float limit,float suma,float* helper,float* restrict release,float harmonic_diff){
    float max_diff = harmonic_diff;
    float* restrict release_ittr = release;
    if(suma<limit){
        for(float* restrict ei = eq;ei<eq+bins;ei++){
            float teq = 1.0;
            float diff = teq - *ei;
            *ei = *ei + diff*(*release);
            release++;
        }
        return;
    }

    float max = 0.0;
    float min = 1.0;

    float* restrict helper_ittr = helper;
    for(float* restrict fft = fft_out;fft<fft_out+bins;fft++,helper_ittr++){
        float ratio = (*fft)/suma;
        if(ratio > max){
            max = ratio;
        }
        if(ratio<min){
            min = ratio;
        }
        *helper_ittr = ratio;

    }
    float diff = max-min;
    //frequency ratio compression
    if(diff > max_diff){
        float ratio = 1.0 - (diff - max_diff);
        float sf = (1.0 - ratio) / bins;
        for(float* restrict r = helper;r<helper+bins;r++){
           *r = (*r) * ratio + sf;
        }

    }



    helper_ittr = helper;
    float* restrict fft_out_ittr = fft_out;
    int* restrict r_ittr = rastoyane;
    for(float* restrict eq_ittr = eq;eq_ittr<eq+bins;eq_ittr++,helper_ittr++,fft_out_ittr++,r_ittr++){
        if(*fft_out_ittr!=0){
            float mult = (*helper_ittr);
            float distance = (*r_ittr);
            float current = *eq_ittr;
            float target = (mult*limit);
            float eq_targ = (target/(*fft_out_ittr));

            if(eq_targ<*eq_ittr || (*release_ittr > 0.4)){
                float diff = eq_targ-current;

                *eq_ittr = *eq_ittr+(diff/distance);
            }else{
                float diff = eq_targ-current;

                *eq_ittr = *eq_ittr+(diff*(*release_ittr));

            }
            release_ittr++;

        }
    }


}
//gain controller
//set to 100 for Pi zero, set to 1024 for normal computers
int agc_lookahead = 300;
struct Gain_Control* gain_control_init(float attack, float release, float target,float noise_th){
    struct Gain_Control* gc = malloc(sizeof(struct Gain_Control));
    gc->attack = attack;
    gc->release = release;
    gc->target = target;
    gc->RMS_sum = malloc(sizeof(float)*agc_lookahead);
    memset(gc->RMS_sum,0,sizeof(float)*agc_lookahead);
    gc->gain = 1.0;
    gc->noise_th = noise_th;
    return gc;
}
void set_gain_control(struct Gain_Control* gc, float target,float attack, float release,float noise_th){
    gc->attack = attack;
    gc->release = release;
    gc->target = target;
    gc->noise_th = noise_th;
}

void gain_control(struct Gain_Control* gc, float* levo, float* pravo){
    float target = gc->target;
    float attack = gc->attack;
    float release = gc->release;
    float left = *levo;
    float right = *pravo;
    float stereo = left-right;
    float sum = ((left+right)+stereo)/2.0;
    sum = sum*sum;

    float run_sum = 0;
    for(float* restrict i = gc->RMS_sum + (agc_lookahead-1);i>gc->RMS_sum;i--){
        float val = *(i-1);
        run_sum = run_sum + val;
        *i = val;
    }
    run_sum = run_sum + sum;
    *(gc->RMS_sum) = sum;




    float rms_val_ps = run_sum/agc_lookahead;

    float rms_val = sqrtf(rms_val_ps);

    if(rms_val<gc->noise_th){
        if(gc->gain<0.9){
            gc->gain = gc->gain + release;
        }
        if(gc->gain>1.1){
            gc->gain = gc->gain - release;
        }
    }else{
        if(rms_val*gc->gain > target){
            gc->gain = gc->gain - attack;
        }else{
            gc->gain = gc->gain + release;
        }
    }


    *levo = left*gc->gain;
    *pravo = right*gc->gain;
}

void free_gain_control(struct Gain_Control *gc){
    free(gc->RMS_sum);
    free(gc);
}
