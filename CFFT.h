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
    static complex<float> * compute(complex<float> * x, int n, int P, int Q, int threads);
private:
    static complex<float> *sliceVector(complex<float>* x, int start, int n, int stride);
    static complex<float> *mergeVectors(complex<float> * a, complex<float> *b, int n);
    static complex<float> twiddle(int N, int k);
    static complex<float> *copyVector(complex<float>* a, int n);
    static void tritconvert(int i, int m, int* trit);
    static void tritrev(int *trit,int m);
    static int trittoint(int *trit,int m);
    static void printVector(std::complex<float> * vector, int n);
    static void deleteMatrix(complex<float>** a, int n);
    static void compute2radix(complex<float> * x, int n);
    static void compute3radix(complex<float> * x, int n);
    static void reorder3Radix(complex<float>*&x, int n);
    static map<pair<int,int>,complex<float> > m_twiddleFactors;
    
};

#endif	/* CFFT_H */

