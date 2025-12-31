FLAGS=-O3 -mfpmath=both -march=native
PI_FLAGS=-mfpu=vfp -O3 -march=armv6zk -mtune=arm1176jzf-s
FFT=./FFT/FFT.c
ALSA=./alsa_pipe/main.c
DEQ=./EQ_max_algorithm/maxim.c

all:
	$(CC) audio.c $(FFT) $(DEQ) $(ALSA) $(FLAGS)  -lm -lasound -Wall -o vtkradio
