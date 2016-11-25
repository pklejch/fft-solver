/* 
 * File:   CFFT.h
 * Author: Petr Klejch (klejcpet@fit.cvut.cz)
 *
 * Created on 12. března 2016, 16:16
 */

#ifndef CFFT_H
#define	CFFT_H

#include <complex>
#include <map>
using namespace std;

class CFFT {
public:
    CFFT();
    static complex<float> * compute(float * reals,float * imags, int n, int P, int Q);
private:
    static complex<float> *mergeVectors(complex<float> * a, complex<float> *b, int n);
    static complex<float> *copyVector(complex<float>* a, int n);
    static void deleteMatrix(complex<float>** a, int n);
    
};

#endif	/* CFFT_H */

