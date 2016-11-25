/* 
 * File:   CFFT.cpp
 * Author: Petr Klejch (klejcpet@fit.cvut.cz)
 * 
 * Created on 12. března 2016, 16:16
 */

#include "CFFT.h"
#include <iostream>
#include <math.h>
#include <omp.h>
#include <stdexcept>
#include <cstdio>

//komplexni cisla pro CUDA
#include <thrust/complex.h>



// CUDA Includes
#include <cuda.h>
#include <cuda_runtime.h>


//limit od kdy se dela DFT pro radix2-FFT
#define LIMIT 32


__global__ void transpose(const float *reals, const float * imags, float *outputReals, float * outputImags, int Q, int P);
__global__ void multiplyTwiddle(float *reals, float * imags, int Q, int P);
__global__ void computeFFT1(float *reals, float * imags, float *tempReal, float* tempImag, float* tempReal2, float* tempImag2,int* trit,float * dftReal, float * dftImag, int Q, int P, bool swapped, int m);
__global__ void computeFFT2(float *reals, float * imags, float *tempReal, float* tempImag, float* tempReal2, float* tempImag2,int* trit,float * dftReal, float * dftImag, int Q, int P, bool swapped, int m);
__device__ void compute2radix(float* reals, float *imags, float *tempReal, float* tempImag, float* tempReal2, float* tempImag2,float * dftReal, float * dftImag, int n, int tid,int size);
__device__ void sliceVector(float* x,float * output, int start, int n, int stride);
__device__ void reorder3Radix(float* reals, float *imags, int * trit, int n, int tid);
__device__ float twiddleReal(int N, int k);
__device__ float twiddleImag(int N, int k);
__device__ thrust::complex<float> twiddle(int N, int k);
__device__ void compute3radix(float* reals, float * imags, int n);
__device__ void computeDFT(float * reals, float* imags, float * tempReal, float * tempImag, int n);
void printVector(float * vector, int n);

__device__ void computeDFT(float * reals, float* imags,float *tempReal, float * tempImag, int n){

    float sumReal = 0;
    float sumImag = 0;
    float angle;
    //iterace pres prvky vstupniho pole
    for (int k=0; k < n; k++){
        sumReal = 0;
        sumImag = 0;

        //iterace pres prvky vstupniho pole
        for (int t=0; t < n;t++){

            angle = 2*M_PI*t*k/n;
            sumReal+= (reals[t]* __cosf(angle)) + (imags[t] * __sinf(angle));
            sumImag+= (-reals[t]* __sinf(angle)) + (imags[t] * __cosf(angle));
        }

        //ulozeni hodnot
        tempReal[k]=sumReal;
        tempImag[k]=sumImag;
    }
}

//P*Q vlaken celkem
//Q- sloupce puvodni matice
//P - radky puvodni matice
//OK - 100%
__global__ void transpose(const float *reals, const float * imags, float *outputReals, float * outputImags, int Q, int P) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;

    if (tid<P*Q) {
        int X = tid % Q;
        int Y = tid / Q;

        //(x,y) => (y,x)
        //printf("thread %d: swapping %d ---> %d \n",tid,tid,(P*X) + Y);
        outputReals[(P*X) + Y] = reals[tid];
        outputImags[(P*X) + Y] = imags[tid];
    }
}

//P*Q vlaken celkem
//Q - sloupce
//P - radky
//OK - 100%
__global__ void multiplyTwiddle(float *reals, float * imags, int Q, int P) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid<P*Q) {
        int X = tid % Q;
        int Y = tid / Q;
        //printf("thread %d mutlipling element %d by [%d,%d] \n",tid,tid,X,Y);
        thrust::complex<float> t = thrust::complex<float>(reals[tid],imags[tid]) * twiddle(P*Q,X*Y);
        reals[tid] = t.real();
        imags[tid] = t.imag();
    }
}


//celkem P vlaken
__global__ void computeFFT1(float *reals, float * imags,float *tempReal, float* tempImag, float* tempReal2, float* tempImag2,int *trit,float * dftReal, float * dftImag, int Q, int P, bool swapped, int m) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid < P){

        if(!swapped){
            //sloupce jsou mocnina 2
        	//printf("thread %d: Running FFT1-radix2, at %d of length %d \n",tid,tid*Q,Q);
            compute2radix(&(reals[tid*Q]),&(imags[tid*Q]),&(tempReal[tid*Q]),&(tempImag[tid*Q]),&(tempReal2[tid*Q]),&(tempImag2[tid*Q]),&(dftReal[tid*Q]),&(dftImag[tid*Q]),Q,tid,Q);
        }else{
            if(Q>3){
                //sloupce jsou mocnina 3
            	//printf("thread %d: Running FFT1-radix3, at %d of length %d \n",tid,tid*Q,Q);
            	//printf("thread %d: trit at %d of len %d\n",tid,tid*m,m);
                compute3radix(&(reals[tid*Q]),&(imags[tid*Q]),Q);
                reorder3Radix(&(reals[tid*Q]),&(imags[tid*Q]),&(trit[tid*m]),Q,tid);
            }
        }
    }
}

//celkem Q vlaken
__global__ void computeFFT2(float *reals, float * imags,float *tempReal, float* tempImag, float* tempReal2, float* tempImag2,int *trit,float * dftReal, float * dftImag, int Q, int P, bool swapped, int m) {
    int tid = threadIdx.x + blockIdx.x * blockDim.x;
    if (tid < Q){

        if (!swapped){
            if(P>3){
            	//printf("thread %d: Running FFT2-radix3, at %d of length %d \n",tid,tid*P,P);
                //sloupce jsou mocnina 3
                compute3radix(&(reals[tid*P]),&(imags[tid*P]),P);
                reorder3Radix(&(reals[tid*P]),&(imags[tid*P]),&(trit[tid*m]),P,tid);
            }
        }else{
            //sloupce jsou mocnina 2
        	//printf("thread %d: Running FFT2-radix2, at %d of length %d \n",tid,tid*P,P);
        	compute2radix(&(reals[tid*P]),&(imags[tid*P]),&(tempReal[tid*P]),&(tempImag[tid*P]),&(tempReal2[tid*P]),&(tempImag2[tid*P]),&(dftReal[tid*P]),&(dftImag[tid*P]),P,tid,P);
        }
    }
}

void checkError( cudaError_t cudaerr){
    if (cudaerr != CUDA_SUCCESS)
        printf("cuda error \"%s\".\n", cudaGetErrorString(cudaerr));
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
void printVector(float * vector, int n){
    for (int i = 0; i < n; i++){
        std::cout<< vector[i] << std::endl; 
    }
}
void printAsMatrix(float * vector, float* vector2, int n, int P, int Q){
	for(int i=0; i<P;i++){
		for(int j=0; j<Q; j++){
			std::cout<<"("<< vector[i*Q+j] << ", "<< vector2[i*Q+j]<<") ";
		}
		std::cout<<endl;
	}
}

void printAsMatrix2(float * vector, float* vector2, int n, int P, int Q){
	for(int i=0; i<Q;i++){
		for(int j=0; j<P; j++){
			std::cout<<"("<< vector[i*P+j] << ", "<< vector2[i*P+j]<<") ";
		}
		std::cout<<endl;
	}
}
void printVectors(float * vector, float * vector2, int n){
    for (int i = 0; i < n; i++){
        std::cout<< "("<<vector[i] << ", "<< vector2[i]<<")"<< std::endl;
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
complex<float> * CFFT::compute(float* reals, float* imags, int n, int P, int Q){
    
	cudaError_t cudaerr;
    //vytvorim novy vystupni vektor
    complex<float>* output2 = new complex<float>[n];
    

    bool swapped = false;
    size_t tempDev_size;
    size_t tritDev_size;
    size_t tempDftDev_size;
    int m;
    //normalizuji Q a P, aby matice byla vzdy na vysku
    if(Q>P){
        int t=Q;
        Q=P;
        P=t;
        swapped = true;
        tempDev_size = sizeof(float) * size_t(Q*P);
        tempDftDev_size = sizeof(float) * LIMIT*Q;
        //pro P vlaken
        m=round(log(Q)/log(3));
        tritDev_size = sizeof(int) * m * P;
    }else{
    	tempDev_size = sizeof(float) * size_t(P*Q);
    	tempDftDev_size = sizeof(float) *size_t(LIMIT*P);
    	// pro Q vlaken
    	m = round(log(P)/log(3));
    	tritDev_size = sizeof(int) * m * Q;
    }

    cout << "Q="<< Q << ", P="<<P<<endl;
    if (swapped){
        cout << "Swapped nastaven na true."<<endl;
    }else{
        cout << "Swapped nastaven na false."<<endl;
    }
    ///////////////////////////////////////////////////////////////////
    //alokace pomocne matice pro transpozici - P sloupcu, Q radku
    float* transposeDevReal;
    float* transposeDevImag;
    const size_t transposeDev_size = sizeof(float) * size_t(P*Q);
    
    //cout << "Alokace mista pro transponovanou matici "<<endl;
    cudaerr= cudaMalloc((void **)&transposeDevReal, transposeDev_size);
    checkError(cudaerr);
    cudaerr= cudaMalloc((void **)&transposeDevImag, transposeDev_size);
    checkError(cudaerr);
    /////////////////////////////////////////////////////////////////
    
    //pomocne vektory
    float * tempDevReal;
    float * tempDevImag;
    cudaerr= cudaMalloc((void **)&tempDevReal, tempDev_size);
    checkError(cudaerr);
    cudaerr= cudaMalloc((void **)&tempDevImag, tempDev_size);
    checkError(cudaerr);

    float * tempDevReal2;
    float * tempDevImag2;
    cudaerr= cudaMalloc((void **)&tempDevReal2, tempDev_size);
    checkError(cudaerr);
    cudaerr= cudaMalloc((void **)&tempDevImag2, tempDev_size);
    checkError(cudaerr);

    int * tritDev;

    //cout << "Alokace pomocneho pole trit o velikosti "<<tritDev_size/sizeof(int)<<endl;
    cudaerr= cudaMalloc((void **)&tritDev, tritDev_size);
    checkError(cudaerr);



    //cout << "Alokace pomocneho pole o veliksoti "<<tempDev_size/sizeof(float)<<endl;

    float * tempDevDftReal;
    float * tempDevDftImag;
    cudaerr= cudaMalloc((void **)&tempDevDftReal,tempDftDev_size );
    cudaerr= cudaMalloc((void **)&tempDevDftImag,tempDftDev_size );
    checkError(cudaerr);

    //v Q je ulozen mensi faktor
    //v P je ulozen vetsi faktor
    //vytvorim 2D matici o velikosti Q sloupcu, P radku

    //alokace 2D matice
        
    //matice pro pocitani
    float * outputDevReal;
    float * outputDevImag;
    
    //cout << "Alokace matice pro pocitani" << endl;
    const size_t outputDev_size = sizeof(float) * size_t(P*Q);
    cudaerr= cudaMalloc((void **)&outputDevReal, outputDev_size);
    checkError(cudaerr);
    cudaerr= cudaMalloc((void **)&outputDevImag, outputDev_size);
    checkError(cudaerr);

    ///////////////////////////////////////////////////////////////////
    //naplnim matici P*Q
    cout << "Naplneni matice pro pocitani"<<endl;

    //reals - > outputDevReal
    cudaerr= cudaMemcpy(outputDevReal, reals, outputDev_size, cudaMemcpyHostToDevice);
    checkError(cudaerr);
    //imags -> outputDevImag
    cudaerr= cudaMemcpy(outputDevImag, imags, outputDev_size, cudaMemcpyHostToDevice);
    checkError(cudaerr);


    int threadsPerBlock;
    int blocksPerGrid;


    int sizeOfBlock=32;

    cudaEvent_t start, stop;
    float elapsedTime;
    cudaEventCreate( &start ) ;
    cudaEventCreate( &stop ) ;

    cudaDeviceSetLimit(cudaLimitStackSize,25*1024);
    cudaDeviceSetLimit(cudaLimitMallocHeapSize,200000*1000);

    //////////////////////////////
    //KONEC ALOKACE
    /////////////////////////////


    ///////////////////////////////////////////////////////////////////
    //transpozice Output - > transpose


    //P*Q vlaken celkem
    if (n<sizeOfBlock){
        threadsPerBlock = n;
        blocksPerGrid   = 1;
    } else {
        threadsPerBlock = sizeOfBlock;
        blocksPerGrid   = ceil(double(n)/double(threadsPerBlock));
    }
    //cout << "transpozice, pocet bloku: "<<blocksPerGrid<<", pocet vlaken v bloku: "<<threadsPerBlock << endl;
    transpose<<<blocksPerGrid,threadsPerBlock>>>(outputDevReal,outputDevImag,transposeDevReal,transposeDevImag,Q,P);
    cudaThreadSynchronize();



    cudaEventRecord(start,0);
    ///////////////////////////////////////////////////////////////////
    ///KERNEL
    ///////////////////////////////////////////////////////////////////
    //udelam Qkrat FFT o velikosti P (po radcich)

    //Q vlaken
    if (Q<sizeOfBlock){
        threadsPerBlock = Q;
        blocksPerGrid   = 1;
    } else {
        threadsPerBlock = sizeOfBlock;
        blocksPerGrid   = ceil(double(Q)/double(threadsPerBlock));
    }
    cout << "prvni cast FFT, pocet bloku: "<<blocksPerGrid<<", pocet vlaken v bloku: "<<threadsPerBlock << endl;
    computeFFT2<<<blocksPerGrid,threadsPerBlock>>>(transposeDevReal,transposeDevImag,tempDevReal,tempDevImag,tempDevReal2,tempDevImag2,tritDev,tempDevDftReal,tempDevDftImag,Q,P,swapped,m);
    cudaThreadSynchronize();



    ///////////////////////////////////////////////////////////////////
    //KERNEL
    ///////////////////////////////////////////////////////////////////
    //transpozice Trans->Output


    //P*Q vlaken celkem
    if (n<sizeOfBlock){
        threadsPerBlock = n;
        blocksPerGrid   = 1;
    } else {
        threadsPerBlock = sizeOfBlock;
        blocksPerGrid   = ceil(double(n)/double(threadsPerBlock));
    }
    cout << "transpozice, pocet bloku: "<<blocksPerGrid<<", pocet vlaken v bloku: "<<threadsPerBlock << endl;
    transpose<<<blocksPerGrid,threadsPerBlock>>>(transposeDevReal,transposeDevImag,outputDevReal,outputDevImag,P,Q);

    cudaThreadSynchronize();

    ///////////////////////////////////////////////////////////////////
    // KERNEL
    /////////////////////////////////////////////////////////////////
    //pronasobeni twiddle faktory
    //P*Q vlaken
    if (n<512){
        threadsPerBlock = n;
        blocksPerGrid   = 1;
    } else {
        threadsPerBlock = sizeOfBlock;
        blocksPerGrid   = ceil(double(n)/double(threadsPerBlock));
    }
    cout << "Nasobeni twiddle faktory, pocet bloku: "<<blocksPerGrid<<", pocet vlaken v bloku: "<<threadsPerBlock << endl;
    multiplyTwiddle<<<blocksPerGrid,threadsPerBlock>>>(outputDevReal,outputDevImag,Q,P);

    cudaThreadSynchronize();


    ///////////////////////////////////////////////////////////////////
    // KERNEL
    ///////////////////////////////////////////////////////////////////
    //udelam Pkrat FFT o velikosti Q (po radcich)
    
    //celkem P vlaken
    if (P<sizeOfBlock){
        threadsPerBlock = P;
        blocksPerGrid   = 1;
    } else {
        threadsPerBlock = sizeOfBlock;
        blocksPerGrid   = ceil(double(P)/double(threadsPerBlock));
    }
    cout << "druha cast FFT, pocet bloku: "<<blocksPerGrid<<", pocet vlaken v bloku: "<<threadsPerBlock << endl;
    computeFFT1<<<blocksPerGrid,threadsPerBlock>>>(outputDevReal,outputDevImag,tempDevReal,tempDevImag,tempDevReal2,tempDevImag2,tritDev,tempDevDftReal,tempDevDftImag,Q,P,swapped,m);

    cudaThreadSynchronize();
    cudaEventRecord(stop,0);


    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&elapsedTime,start,stop);
    cout << "Vypocet trval " <<elapsedTime<<" ms"<<endl;



    ///////////////////////////////////////////////////////////////////
    // CUDAMEMCPY
    ///////////////////////////////////////////////////////////////////
    //zapisu vyslednou matici do vektoru po radcich

    

    cudaerr= cudaMemcpy(reals, outputDevReal, transposeDev_size, cudaMemcpyDeviceToHost);

    cudaerr= cudaMemcpy(imags, outputDevImag, transposeDev_size, cudaMemcpyDeviceToHost);

    ////////////////////////////////////////////////////////////////////
    //konec vypoctu
    cout << "Ukladam vysledny vektor"<<endl;


    int k=0;
    for(int j=0;j<Q;j++){ //j <- 0..Q
    	for(int i=0;i<P;i++){ //i <- 0..P
    		output2[k] = complex<float>(reals[i*Q+j],imags[i*Q+j]);
    		k++;
    	}
    }

    //smazu nepotrebne vektory
    cudaFree(transposeDevReal);
    cudaFree(transposeDevImag);
    cudaFree(tritDev);
    cudaFree(tempDevReal);
    cudaFree(tempDevImag);
    cudaFree(tempDevReal2);
    cudaFree(tempDevImag2);
    cudaFree(outputDevReal);
    cudaFree(outputDevImag);


    return output2;
}

/*
 Spocita FFT pro vektor o velikosti n = 2^k
 */
__device__ void compute2radix(float* reals, float *imags, float *tempReal, float* tempImag, float* tempReal2, float* tempImag2,float *dftReal, float *dftImag, int n,int tid,int size){

	int newLen = n/2;

    //od nejake meze delat DFT
    if(n<=LIMIT){


        computeDFT(reals,imags,dftReal,dftImag,n);
        for(int i=0; i<n; i++){
        	//printf("%f , %f",reals[i],imags[i]);
            reals[i]=dftReal[i];
            imags[i]=dftImag[i];
        }
    }else{
    	//printf("thread %d: newLen = %d\n",tid,newLen);
    	if(newLen<=0){
    		printf("tohle se nesmi stat\n");
    		return;
    	}
    	//printf("thread %d: n = %d\n",tid,n);
        //rozdel


        sliceVector(reals,tempReal,0,n,2); //even real
        sliceVector(imags,tempImag,0,n,2); //even imag
        
        sliceVector(reals,tempReal2,1,n,2); //odd real
        sliceVector(imags,tempImag2,1,n,2); //odd imag
        /*
        if(!tid){
        	printf("after slice: \n");
        	printf("even:\n");
        }
        for(int i=0; i<n/2; i++){
        	printf("thread %d: %f , %f \n",tid,tempReal[i],tempImag[i]);
        }
        if(!tid)
        	printf("odd:\n");
        for(int i=0; i<n/2; i++){
        	printf("thread %d: %f , %f \n",tid,tempReal2[i],tempImag2[i]);
        }*/
        //panuj
        //printf("thread %d:start of temp array at index %d of length %d \n",tid,newLen,newLen);
        //int len = &(tempReal[newLen]) - tempReal;
        //printf("thread %d: lenght of array %d \n",tid,len);

        compute2radix(tempReal, tempImag, &(tempReal[newLen]), &(tempImag[newLen]),&(tempReal2[newLen]), &(tempImag2[newLen]),dftReal,dftImag,newLen,tid,size);
        compute2radix(tempReal2,tempImag2,&(tempReal[newLen]), &(tempImag[newLen]),&(tempReal2[newLen]), &(tempImag2[newLen]),dftReal,dftImag,newLen,tid,size);
        

        //sloz dohromady
        for(int k=0; k < (newLen); k++ ){
            //twiddle factor
            
            thrust::complex<float> t = twiddle(n,k) * thrust::complex<float>(tempReal2[k],tempImag2[k]);

            thrust::complex<float> t2 = thrust::complex<float>(tempReal[k],tempImag[k]) + t;

            //printf("start of array %p \n",&(reals[k]));
            reals[k] = t2.real();
            imags[k] = t2.imag();
            
            thrust::complex<float> t3 = thrust::complex<float>(tempReal[k],tempImag[k]) - t;

            reals[k+newLen] = t3.real();
            imags[k+newLen] = t3.imag();


        }
    }
}
__device__ float twiddleReal(int N, int k){
    return __cosf(k*2*M_PI/(float)N);
}

__device__ float twiddleImag(int N, int k){
    return -__sinf(k*2*M_PI/(float)N);
}

__device__ thrust::complex<float> twiddle(int N, int k){
	thrust::complex<float> output(__cosf(k*2*M_PI/(float)N),-__sinf(k*2*M_PI/(float)N));
    return output;
}
/*
 Spocita FFT o velikosti n=3^k
 */
__device__ void compute3radix(float* reals, float * imags, int n){
    int k1, N1, N2;
    thrust::complex<float> butterfly[3];
    N1=3;
    N2=n/3;
    for (int n2=0;n2<N2;n2++){

    	butterfly[0] = (thrust::complex<float>(reals[n2],imags[n2]) + thrust::complex<float> (reals[N2 + n2],imags[N2 + n2]) + thrust::complex<float> (reals[2*N2 + n2],imags[2*N2 + n2]));
        butterfly[1] = twiddle(n,n2)*(thrust::complex<float>(reals[n2],imags[n2]) + (thrust::complex<float>(reals[N2 + n2],imags[N2 + n2]) * twiddle(3,1)) + (thrust::complex<float>(reals[2*N2 + n2],imags[2*N2 + n2])*twiddle(3,2)));
        butterfly[2] = twiddle(n,(2*n2)) * (thrust::complex<float>(reals[n2],imags[n2])+(thrust::complex<float>(reals[N2+n2],imags[N2+n2])*twiddle(3,2)) + (thrust::complex<float>(reals[2*N2 + n2],imags[2*N2 + n2])*twiddle(3,1)));



        for(k1=0;k1<N1;k1++){
            reals[n2+N2*k1]=butterfly[k1].real();
            imags[n2+N2*k1]=butterfly[k1].imag();
        }
        
    }
    
    if (N2!=1){
        for(k1=0;k1<N1;k1++){
            compute3radix(&(reals[N2*k1]),&(imags[N2*k1]), N2);
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


__device__ void tritconvert(int i, int m, int* trit)
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

__device__ void tritrev(int *trit,int m)	//used in reordering loop (rev the trits)
{
	int temp,i,j;
	for(i=0,j=m-1;i<m/2;i++,j--)
	{
		temp = trit[i];
		trit[i] = trit[j];
		trit[j] = temp;
	}
}


__device__ int trittoint(int *trit,int m)	//used in reordering loop(convert to int)
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
__device__ void reorder3Radix(float* reals, float *imags,int *trit, int n, int tid) {
    
   
    //exponent
    int m=round(__logf(n)/__logf(3));

    int j;
    
    float tempReal;
    float tempImag;
    
    for(int i=0;i<n;i++)
    {
            tritconvert(i,m,trit);		//convert to trits
            tritrev(trit,m);			//reverse the trit
            j = trittoint(trit,m);		//convert to integer

            if(j>i)						//swap
            {
                tempReal = reals[i];
                tempImag = imags[i];
                
                reals[i]=reals[j];
                imags[i]=imags[j];
                
                reals[j]=tempReal;
                imags[j]=tempImag;                        
            }
    }
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

__device__ void sliceVector(float* x,float * output, int start, int n, int stride){
    //int j;
    int j=0;
    for(int i=start; i < n; i+=stride){
        output[j]=x[i];
        j++;
    }
}
