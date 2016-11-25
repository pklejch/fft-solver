CC=g++ -Wall -pedantic -fopenmp -O3
ICC=icc -mmic -openmp -Ofast
SOURCES=$(wildcard *.cpp)
OBJECTS=$(wildcard *.o)

all: pap_fft

pap_fft_vect:
	$(CC) -ftree-vectorize -ffast-math -fopt-info-vec -ftree-vectorizer-verbose=6 -o $@ $(SOURCES)
pap_fft:
	$(CC) -o $@ $(SOURCES)
pap_fft_phi:
	$(ICC) -vec-report5 -o $@ $(SOURCES)

pap_fft_phi_novect:
	$(ICC) -no-vec -o $@ $(SOURCES)

cuda_star:
	nvcc -O3 -gencode arch=compute_30,code=sm_30 --use_fast_math -g -G -c -o CFFT.o CFFT.cu
	nvcc -O3 -gencode arch=compute_30,code=sm_30 --use_fast_math -g -G -c -o CParser.o CParser.cpp
	nvcc -O3 -gencode arch=compute_30,code=sm_30 --use_fast_math -g -G -c -o main.o main.cpp
	nvcc -O3 -gencode arch=compute_30,code=sm_30 --use_fast_math -g -G --ptxas-options=-v -o $@ $(OBJECTS)

clean:
	rm pap_fft pap_fft_phi
