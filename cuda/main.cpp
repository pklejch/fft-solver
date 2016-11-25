/* 
 * File:   main.cpp
 * Author: petr
 *
 * Created on 12. března 2016, 11:12
 */

#include <cstdlib>
#include <complex>
#include <iostream>
#include <fstream>
#include "CFFT.h"
#include "CParser.h"
#include <string>



using namespace std;

void printVector(std::complex<float> * vector, int n){
    for (int i = 0; i < n; i++){
        std::cout<< vector[i] << std::endl; 
    }
}


void writeVector(std::complex<float> * vector, int n, string file){
    ofstream output(file.c_str());
    for (int i = 0; i < n; i++){
        output<< vector[i] << std::endl; 
    }
}

void findFactors(int N, int& P, int & Q){
    int k=0;
    int originN = N;
    while(N%3 == 0){
        k++;
        N=N/3;
    }
    P=(int)pow(3,k);
    Q=originN/P;
    int testQ = Q;
    
    while(true){
        if(testQ%2 == 0){
            testQ=testQ/2;
        }else{
            if(testQ == 1){
                break;
            }else{
                cout << "Velikost vstupniho vektoru neni ve tvaru n=2^k1*3^k2" <<endl;
                exit(EXIT_FAILURE);
            }
        }
    }
}
/*
 * 
 */
int main(int argc, char** argv) {
    int N=0,P=0,Q=0;
  
    float * reals;
    float * imags;
    cout << "Nacitani vstupniho vektoru..."<<endl;
    CParser::parse(reals,imags,N);
  
    findFactors(N,P,Q);
  
    cout << "Velikost vektoru je " << N << ", velikost odpovidajici matice je "<< Q << "x"<<P<<endl;
    cout << "Startuji vypocet..."<<endl;
    complex<float> * output2 = CFFT::compute(reals,imags,N,P,Q);
    
    cout << "Zapis vysledneho vektoru..."<<endl;
    writeVector(output2,N,"FFT_out.txt");
    
    delete [] reals;
    delete [] imags;
    delete [] output2;
    return 0;
}

