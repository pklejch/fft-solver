/* 
 * File:   CFFT.cpp
 * Author: Petr Klejch (klejcpet@fit.cvut.cz)
 * 
 * Created on 12. března 2016, 16:16
 */

#include "CFFT.h"
#include "CDFT.h"
#include <iostream>
#include <math.h>
#include <omp.h>
#include <stdexcept>

map<pair<int,int>,complex<float> > CFFT::m_twiddleFactors;

CFFT::CFFT() {
}


void printMatrix(complex<float>** a, int rows,int cols){
    for(int i=0; i<rows; i++){
        for(int j=0; j<cols; j++){
            cout << a[j][i] << " ";
        }
        cout << endl;
    }
    cout << endl;
}
void CFFT::printVector(std::complex<float> * vector, int n){
    for (int i = 0; i < n; i++){
        std::cout<< vector[i] << std::endl; 
    }
}

void CFFT::deleteMatrix(complex<float>** a, int n){
    for(int i=0; i<n; i++){
        delete [] a[i];
    }
    delete [] a;
}
/*
 Spocita FFT pro vektor velikosti n = 2^k1*3^k2
 */
complex<float> * CFFT::compute(complex<float>* x, int n, int P, int Q, int threads){
    

    
    //matice pro pocitani
    complex<float> **output;
    
    //vytvorim novy vystupni vektor
    complex<float>* output2 = new complex<float>[n];
    
    //vytvoreni transponovane matice P sloupcu, Q radku
    //alokace

    
    bool swapped = false;
    //normalizuji Q a P, aby matice byla vzdy na vysku
    if(Q>P){
        int t=Q;
        Q=P;
        P=t;
        swapped = true;
    }
    
    complex<float>** __restrict__ transpose = new complex<float>*[P];
    //nactu startovni cas
    double start = omp_get_wtime();
    omp_set_num_threads(threads);
    #pragma omp parallel shared(output,output2,transpose)
    {
        
        #pragma omp for
        for(int i=0; i<P; i++){
            transpose[i]=new complex<float>[Q];
        }
        //v Q je ulozen mensi faktor
        //v P je ulozen vetsi faktor
        //vytvorim 2D matici o velikosti Q sloupcu, P radku

        //alokace 2D matice
        output = new complex<float>*[Q];
        #pragma omp for
        for(int i=0;i<Q;i++){
            output[i] = new complex<float>[P];
        }

        //naplnim matici P*Q
        # pragma omp for collapse (2)
        for(int i=0;i<Q;i++){
            for(int j=0;j<P;j++){
                output[i][j]=x[i+Q*j];
            }
        }
        //udelam Qkrat FFT o velikosti P (po sloupcich)
        if(swapped){
            #pragma omp for
            for(int i=0;i<Q;i++){
                //sloupce jsou mocnina 2
                CFFT::compute2radix(output[i],P);
            }
        }else{
            if(P>3){
                #pragma omp for
                for(int i=0;i<Q;i++){
                    //sloupce jsou mocnina 3
                    CFFT::compute3radix(output[i],P);
                    CFFT::reorder3Radix(output[i],P);
                }
            }
        }

        //pronasobeni twiddle faktory
        # pragma omp for collapse (2)
        for(int i=0;i<Q;i++){
            for(int j=0;j<P;j++){
                output[i][j] *= CFFT::twiddle(n,i*j);
            }
        }
        //transpozice

        //naplneni transponovane matice
        # pragma omp for collapse (2)
        for(int i=0;i<Q;i++){
            for(int j=0;j<P;j++){
                transpose[j][i]=output[i][j];
            }
        }

        //udelam Pkrat FFT o velikosti Q (po sloupcich)
        
        if (swapped){
            if(Q>3){
                #pragma omp for
                for(int i=0;i<P;i++){
                    //sloupce jsou mocnina 3
                    CFFT::compute3radix(transpose[i],Q);
                    CFFT::reorder3Radix(transpose[i],Q);
                }
            }
        }else{
            #pragma omp for
            for(int i=0;i<P;i++){
                //sloupce jsou mocnina 2
                CFFT::compute2radix(transpose[i],Q);
            }
        }
        
        //zapisu vyslednou matici do vektoru po radcich
        # pragma omp for collapse (2)
        for(int i=0;i<Q;i++){
            for(int j=0;j<P;j++){
                output2[P*i+j]=transpose[j][i];
            }
        }
    }
    //konec vypoctu
    double end = omp_get_wtime();
    
    //smazu nepotrebne matice
    CFFT::deleteMatrix(output,Q);
    CFFT::deleteMatrix(transpose,P);
    cout << "Vypocet trval "<<end-start<<"s."<<endl;       
    return output2;
}

/*
 Spocita FFT pro vektor o velikosti n = 2^k
 */
void CFFT::compute2radix(complex<float>* x, int n){

    //TODO od nejake meze delat DFT
    if(n<=2){
        complex<float> * temp;
        temp=CDFT::compute(x,n);
        for(int i=0; i<n; i++){
            x[i]=temp[i];
        }
        delete [] temp;
    }else{
        complex<float>* even;
        complex<float>* odd;
        
        //rozdel
        even = sliceVector(x,0,n,2);
        odd = sliceVector(x,1,n,2);
        
        
        //panuj
        compute2radix(even,n/2);
        compute2radix(odd,n/2);
        

        //sloz dohromady
        for(int k=0; k < (n/2); k++ ){
            //twiddle factor
            
            complex<float> t = CFFT::twiddle(n,k) * odd[k];
            x[k] = even[k] + t;
            x[k+n/2] = even[k] - t;
        }
        delete [] even;
        delete [] odd;
    }
}
complex<float> CFFT::twiddle(int N, int k){
    complex<float> output;
    output.real()=cos(k*2*M_PI/(double)N);
    output.imag()=-sin(k*2*M_PI/(double)N);
    return output;
}
/*
 Spocita FFT o velikosti n=3^k
 */
void CFFT::compute3radix(complex<float>* x, int n){
    int k1, N1, N2;
    complex<float> butterfly[3];
    N1=3;
    N2=n/3;
    for (int n2=0;n2<N2;n2++){
        butterfly[0] = (x[n2] + x[N2 + n2] + x[2*N2 + n2]);
        
        butterfly[1] = CFFT::twiddle(n, n2)*(x[n2] + (x[N2 + n2]*CFFT::twiddle(3,1)) + (x[2*N2 + n2]*CFFT::twiddle(3,2)));
        
        butterfly[2] = CFFT::twiddle(n, (2*n2))*(x[n2] + (x[N2 + n2]*CFFT::twiddle(3,2)) + (x[2*N2 + n2]*CFFT::twiddle(3,1)));

            
        for(k1=0;k1<N1;k1++){
            x[n2+N2*k1]=butterfly[k1];
        }
        
    }
    
    //TODO - od nejake meze delat DFT
    if (N2!=1){
        for(k1=0;k1<N1;k1++){
            CFFT::compute3radix(&x[N2*k1], N2);
        }
    }
    
}


complex<float>* CFFT::copyVector(complex<float>* a, int n){
    complex<float>* copy = new complex<float>[n];
    for (int i=0; i<n; i++){
        copy[i].real()=a[i].real();
        copy[i].imag()=a[i].imag();
    }
    return copy;
}


void CFFT::tritconvert(int i, int m, int* trit)
//used inreordering loop (conv to trit)
{
	int j=0;
	for(j=0;j<m;j++)
		trit[j]=0;
	
	j=0;
	while(i>0)
	{
		trit[j]=i%3;
		i = i/3;
		j++;
	}
}

void CFFT::tritrev(int *trit,int m)	//used in reordering loop (rev the trits)
{
	int temp,i,j;
	for(i=0,j=m-1;i<m/2;i++,j--)
	{
		temp = trit[i];
		trit[i] = trit[j];
		trit[j] = temp;
	}
}


int CFFT::trittoint(int *trit,int m)	//used in reordering loop(convert to int)
{
	int i,j=0,z=1;
	for(i=0;i<m;i++)
	{
		j = j + trit[i]*z;
		z*= 3;
	}
	return (j);
}

/*
 Preusporada vystupni vektor do spravneho formatu
 */
void CFFT::reorder3Radix(complex<float>* &x, int n) {
    
   
    //exponent
    int m=round(log(n)/log(3));

    
    int* trit = new int[m]; 
    

    int j;
    
    complex<float> temp;
    for(int i=0;i<n;i++)
    {
            tritconvert(i,m,trit);		//convert to trits
            tritrev(trit,m);			//reverse the trit
            j = trittoint(trit,m);		//convert to integer

            if(j>i)						//swap
            {
                temp = x[i];
                x[i]=x[j];
                x[j]=temp;
                        
            }
    }
    
    delete [] trit;
}



complex<float>* CFFT::mergeVectors(complex<float>* a, complex<float>* b, int n){
    complex<float> * output;
    output = new complex<float>[2*n];
    
    for(int i=0; i < n; i++){
        output[i]=a[i];
        output[i+n]=b[i];
    }
    
    return output;
}

complex<float>* CFFT::sliceVector(complex<float>* x, int start, int n, int stride){
    complex<float> * output;
    output = new complex<float>[n/stride];
    //int j;
    int j=0;
    for(int i=start; i < n; i+=stride){
        output[j]=x[i];
        j++;
    }
    
    return output;
}