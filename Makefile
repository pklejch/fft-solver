CC=g++ -Wall -pedantic -fopenmp -O3
ICC=icc -mmic -openmp -Ofast

SOURCES=$(wildcard *.cpp)

all: pap_fft

pap_fft_vect:
	$(CC) -ftree-vectorize -ffast-math -fopt-info-vec -ftree-vectorizer-verbose=6 -o $@ $(SOURCES)
pap_fft:
	$(CC) -o $@ $(SOURCES)
pap_fft_phi:
	$(ICC) -vec-report5 -o $@ $(SOURCES)

pap_fft_phi_novect:
	$(ICC) -no-vec -o $@ $(SOURCES)

clean:
	rm pap_fft pap_fft_phi
