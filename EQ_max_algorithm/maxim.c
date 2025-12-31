#include "maxim.h"
#include <stdlib.h>
#include <math.h>


double find_amp(double *fft_out, int bins){
    double suma = 0;
    for(double* restrict ittr = fft_out; ittr < fft_out + bins; ittr++){
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
void adjust_eq(double *eq,double* fft_out,int* rastoyane, int bins, double limit,double suma,double* helper,double release,double harmonic_diff){
    double max_diff = harmonic_diff;

    if(suma<limit){
        for(double* restrict ei = eq;ei<eq+bins;ei++){
            double teq = 1.0;
            double diff = teq - *ei;
            *ei = *ei + diff*release;
        }
        return;
    }

    double max = 0.0;
    double min = 1.0;

    double* restrict helper_ittr = helper;
    for(double* restrict fft = fft_out;fft<fft_out+bins;fft++,helper_ittr++){
        double ratio = (*fft)/suma;
        if(ratio > max){
            max = ratio;
        }
        if(ratio<min){
            min = ratio;
        }
        *helper_ittr = ratio;

    }
    double diff = max-min;
    //frequency ratio compression
    if(diff > max_diff){
        double ratio = 1.0 - (diff - max_diff);
        double sf = (1.0 - ratio) / bins;
        for(double* restrict r = helper;r<helper+bins;r++){
           *r = (*r) * ratio + sf;
        }

    }



    helper_ittr = helper;
    double* restrict fft_out_ittr = fft_out;
    int* restrict r_ittr = rastoyane;
    for(double* restrict eq_ittr = eq;eq_ittr<eq+bins;eq_ittr++,helper_ittr++,fft_out_ittr++,r_ittr++){
        if(*fft_out_ittr!=0){
            double mult = (*helper_ittr);
            double target = (mult*limit);
            double eq_targ = (target/(*fft_out_ittr));
            double distance = (*r_ittr);
            double current = *eq_ittr;

            double diff = eq_targ-current;

            *eq_ittr = *eq_ittr+(diff/distance);

        }
    }


}
//gain controller

struct Gain_Control* gain_control_init(double attack, double release, double target,double noise_th){
    struct Gain_Control* gc = malloc(sizeof(struct Gain_Control));
    gc->attack = attack;
    gc->release = release;
    gc->target = target;
    gc->RMS_sum = malloc(sizeof(double)*1024);
    gc->gain = 1.0;
    gc->noise_th = noise_th;
    return gc;
}
void set_gain_control(struct Gain_Control* gc, double target,double attack, double release,double noise_th){
    gc->attack = attack;
    gc->release = release;
    gc->target = target;
    gc->noise_th = noise_th;
}

void gain_control(struct Gain_Control* gc, double* levo, double* pravo){
    double target = gc->target;
    double attack = gc->attack;
    double release = gc->release;
    double left = *levo;
    double right = *pravo;
    double sum = (left+right);
    sum = sum*sum;

    double run_sum = 0;
    for(int i = 1023;i>0;i--){
        double val = gc->RMS_sum[i-1];
        run_sum = run_sum + val;
        gc->RMS_sum[i] = val;
    }
    run_sum = run_sum + sum;
    gc->RMS_sum[0] = sum;




    double rms_val_ps = run_sum/1024.0;

    double rms_val = sqrt(rms_val_ps);

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
