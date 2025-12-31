FLAGS=-O3 -mfpmath=both -march=native
PI_FLAGS=-O3 -march=native -mfpu=vfp -mfloat-abi=hard
OPT_FLAGS=-ffast-math -funroll-loops -ftree-vectorize -flto -fno-signed-zeros -fno-trapping-math
FFT=./FFT/FFT.c
ALSA=./alsa_pipe/main.c
DEQ=./EQ_max_algorithm/maxim.c



all:
	$(CC) audio.c $(FFT) $(DEQ) $(ALSA) $(FLAGS) $(OPT_FLAGS) -g -lm -lasound -Wall -o vtkradio
