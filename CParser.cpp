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

complex<float> * CParser::parse(int &N){
    complex<float> * output;
    int len;
    ifstream infile("input.txt");
    double real,imag;
    
    infile>>len;
    output = new complex<float>[len];
    int i = 0;
    while (infile >> real >> imag)
    {
        output[i] = complex<float>(real, imag);
        i++;
    }
    
    N=len;
    return output;
    
}
