FLAGS=-O3 -mfpmath=both -march=native
FFT=./FFT/FFT.c
ALSA=./alsa_pipe/main.c
DEQ=./EQ_max_algorithm/maxim.c

all:
	$(CC) audio.c $(FFT) $(DEQ) $(ALSA) $(FLAGS)  -lm -lasound -Wall -o vtkradio
