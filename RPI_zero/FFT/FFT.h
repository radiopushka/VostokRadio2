#ifndef FFT_RESAMPLE
#define FFT_RESAMPLE

struct FFT_rsmp{
    //low pass filter for each bin
    float alpha;
    float nalpha;

    float hpf_alpha;
    float hpf_nalpha;

    int bins;

    //high pass filter to remove DC
    float pre_high_pass;

    //ring buffer

    int ring_buffer_size;
    int* count_buff;
    int* count_buff_i;
    int* count_buff_e;
    int* rastoyanee;

    float* lpf_r;
    float** lpf_r_v;
    float** lpf_r_i;
    float* lpf_i;
    float** lpf_i_v;
    float** lpf_i_i;

    float** lpf_i_l;
    float** lpf_r_l;


    float* dsp_array;

    int length;
    int counter;

    //amplitude
    float* amplitude;
    float** amplitude_ring;
    float** amplitude_ring_end;
    float** amplitude_ring_i;

    int OFbins;
    float OF_scale;//minimum OF scaling for minimal distortion
};

struct FFT_rsmp *FFT_resample_init(int bins,int OF_bins,int ring_buffer_delay, float fs, float fend,float srate);
int get_resamp_size(struct FFT_rsmp *rsmp);
float* resamp_pre_process(struct FFT_rsmp *rsmp, float in,float* eq);

// get the audio signal after fourier resampling
float resamp_get_signal(struct FFT_rsmp *rsmp, float* eq);

void free_resamp(struct FFT_rsmp *rsmp);

#endif
