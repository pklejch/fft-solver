/* 
 * File:   CDFT.h
 * Author: petr
 *
 * Created on 12. března 2016, 15:13
 */

#include <complex>

#ifndef CDFT_H
#define	CDFT_H

using namespace std;
/*
 Staticka trida pro vypocet DFT 
 */
class CDFT {
public:
    CDFT();
    static complex<float> *compute(complex<float> * x, int n);
};

#endif	/* CDFT_H */

