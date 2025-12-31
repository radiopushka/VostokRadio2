#ifndef MAXIM_H
#define MAXIM_H

//FFT resampling based limiter and bandwidth expander
double find_amp(double* fft_out,int bins);
void adjust_eq(double* eq,double* fft_out,int* rastoyane,int bins,double limit,double suma,double* helper,double release,double harmonic_diff);

//gain controller
struct Gain_Control{
    double* RMS_sum;
    double attack;
    double release;
    double target;
    double noise_th;
    double gain;
};
struct Gain_Control* gain_control_init(double attack, double release, double target,double noise_th);
void gain_control(struct Gain_Control* gc, double* levo,double* pravo);
void set_gain_control(struct Gain_Control* gc, double target,double attack, double release,double noise_th);
void free_gain_control(struct Gain_Control* gc);

#endif
