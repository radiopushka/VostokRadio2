#ifndef FFT_RESAMPLE
#define FFT_RESAMPLE

struct FFT_rsmp{
    //low pass filter for each bin
    double alpha;
    double nalpha;

    double hpf_alpha;
    double hpf_nalpha;

    int bins;

    //high pass filter to remove DC
    double pre_high_pass;

    //ring buffer

    int ring_buffer_size;
    int* count_buff;
    int* count_buff_i;
    int* count_buff_e;
    int* rastoyanee;

    double* lpf_r;
    double** lpf_r_v;
    double** lpf_r_i;
    double* lpf_i;
    double** lpf_i_v;
    double** lpf_i_i;

    double** lpf_i_l;
    double** lpf_r_l;


    double* dsp_array;

    int length;
    int counter;

    //amplitude
    double* amplitude;
    double** amplitude_ring;
    double** amplitude_ring_end;
    double** amplitude_ring_i;
};

struct FFT_rsmp *FFT_resample_init(int bins,int ring_buffer_delay, float fs, float fend,float srate);
int get_resamp_size(struct FFT_rsmp *rsmp);
double* resamp_pre_process(struct FFT_rsmp *rsmp, double in,double* eq);

// get the audio signal after fourier resampling
double resamp_get_signal(struct FFT_rsmp *rsmp, double* eq);

void free_resamp(struct FFT_rsmp *rsmp);

#endif
