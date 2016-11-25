/* 
 * File:   CDFT.cpp
 * Author: Petr Klejch (klejcpet@fit.cvut.cz)
 * 
 * Created on 12. března 2016, 15:13
 */


#include <math.h>
#include "CDFT.h"

CDFT::CDFT() {
}

/*
 Spocita Diskretni Fourierovu Transformaci v case O(n^2)
 @return dynamicke pole komlexnich cisel
 */

std::complex<float>* CDFT::compute(complex<float> * x, int n){
    
    //vytvoreni vystupnihof pole
    
    complex<float>* output;
    output = new complex<float>[n];
    double sumReal = 0;
    double sumImag = 0;
    double angle;
    //iterace pres prvky vstupniho pole
    //#pragma omp parallel for
    for (int k=0; k < n; k++){
        sumReal = 0;
        sumImag = 0;

        //iterace pres prvky vstupniho pole    
        for (int t=0; t < n;t++){

            angle = 2*M_PI*t*k/n;
            sumReal+= (x[t].real()* cos(angle)) + (x[t].imag() * sin(angle));
            sumImag+= (-x[t].real()* sin(angle)) + (x[t].imag() * cos(angle));
        }

        //ulozeni hodnot
        output[k].real()=sumReal;
        output[k].imag()=sumImag;
    }
    return output;
}
