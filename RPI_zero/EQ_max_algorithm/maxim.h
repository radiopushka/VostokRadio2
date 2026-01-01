#ifndef MAXIM_H
#define MAXIM_H

//FFT resampling based limiter and bandwidth expander
float find_amp(float* fft_out,int bins);
void adjust_eq(float* eq,float* fft_out,int* rastoyane,int bins,float limit,float suma,float* helper,float* release,float harmonic_diff);

//gain controller
struct Gain_Control{
    float* RMS_sum;
    float attack;
    float release;
    float target;
    float noise_th;
    float gain;
};
struct Gain_Control* gain_control_init(float attack, float release, float target,float noise_th);
void gain_control(struct Gain_Control* gc, float* levo,float* pravo);
void set_gain_control(struct Gain_Control* gc, float target,float attack, float release,float noise_th);
void free_gain_control(struct Gain_Control* gc);

#endif
