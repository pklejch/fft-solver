/* 
 * File:   CParser.cpp
 * Author: Petr Klejch (klejcpet@fit.cvut.cz)
 * 
 * Created on 14. března 2016, 9:30
 */

#include "CParser.h"
#include <fstream>
#include <istream>

CParser::CParser() {
}

void CParser::parse(float* &reals, float* &imags, int &N){    
    int len;
    ifstream infile("input.txt");
    double real,imag;
    
    infile>>len;
    reals = new float[len];
    imags = new float[len];
    
    int i = 0;
    while (infile >> real >> imag)
    {
        reals[i] = real;
        imags[i] = imag;
        i++;
    }
    
    N=len;
    return;
    
}
