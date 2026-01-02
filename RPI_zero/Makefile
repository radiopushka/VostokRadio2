PI_FLAGS=-mfpu=vfp -mfloat-abi=hard
OPT_FLAGS=-O3 -march=native -ffast-math -funroll-loops -ftree-vectorize -flto -fno-signed-zeros -fno-trapping-math
FFT=./FFT/FFT.c
CONFIG=./config_file/config.c
ALSA=./alsa_pipe/main.c
DEQ=./EQ_max_algorithm/maxim.c

# lite version optimized and tested to run on the Raspberry Pi Zero W
# CPU usage of 96-100% is normal

all:
	$(CC) audio.c $(FFT) $(DEQ) $(CONFIG) $(ALSA) $(PI_FLAGS) $(OPT_FLAGS)  -lm -lasound -Wall -o vtkradio
