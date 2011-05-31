/*
 *  WHRTF.cu
 *  WHRTF
 *
 *  Created by Diego Gomes on 23/03/11.
 *  Copyright (c) 2011 __MyCompanyName__. All rights reserved.
 *
 */
 
 // Includes
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <sys/time.h>
#include "runtime.h"

// includes, project
#include <cufft.h>

#define DB8_CONST "db8"

// Global variables
int threadsPerBlock = 8;
int blocksPerGrid = 16;

// ################## Internal Functions

// Find delay
double sumSquaresOfArrayElements(float* hrtf, int vecLength);
double sumArrayElements(float* vec0, int vecLength);
float** getSquaresOfArrayElements(float* vec0, float* vec1, int vecLength);
float** sumColumns(float* vec0, float* vec1, int vecLength);
short* findIndexesGreaterThan(float* vec, int vecLength, double minValue);

// Sparse coefficients
double** leFiltros(char* filtro, int* filtroLength);
int max(int* vec, int vecLength);
float*** computeIpAux(float*** Hp, int hlength, float*** Fp, int flength);
float** calculateGp(float** Rp, int rpLength, float*** Fp, int fpLength, float*** Ip, int ipLength);
void dec_poly(double** h, int hlength, float* sist, int ho1dLength, float** G, int* G_size);
double* getPartialVector(double* vec, int numOfElements);

// CUFFT
// Complex data type
typedef float2 Complex; 
static __device__ __host__ inline Complex ComplexAdd(Complex, Complex);
static __device__ __host__ inline Complex ComplexScale(Complex, float);
static __device__ __host__ inline Complex ComplexMul(Complex, Complex);
static __global__ void ComplexPointwiseMulAndScale(Complex*, const Complex*, int, float);

// Filtering functions
void Convolve(const Complex*, int, const Complex*, int, Complex*);

// Padding functions
int PadData(const Complex*, Complex**, int, const Complex*, Complex**, int);

// declaration, forward
//float* convFFT(float* signal, int signalLength, float* filter, int filterLength);
float* convSimple(float* signal, int signalLength, float* filter, int filterLength);


// Memory management
void cleanHostMemory(void* h);
void cleanDeviceMemory(void* d);

// CUDA Error
void checkCUDAError(const char *msg);


// Extern functions
extern "C" short* findDelay(float** hrtf, int length);
extern "C" float* shiftInParallel(float* vec, int vecLength, short delay, int maxLength);
extern "C" void coef_spars(char* filtro[], int filtroLength, float* ho1d, int ho1dLength, float** G_aux, int* G_size);
extern "C" float* resp_imp(char* filtros[], int numFiltros, double** G, int* G_size, int* resultLength);
extern "C" float* convFFT(float* signal, int signalLength, float* filter, int filterLength);


// ################## Device code

// Find delay
__global__ void VecMult(const float* vec0, const float* vec1, float* result0, float* result1, int vecLength) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < vecLength) {
        result0[i] = (vec0[i] * vec0[i]);
		result1[i] = (vec1[i] * vec1[i]);
	}
}

__global__ void sumColumnsInParallel(float* vec, float* result, int vecLength) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
	if (i < vecLength) {
		float sum = 0.0;
		for (int j = 0; j < i; j++) {
			sum += vec[j];
		}		
		result[i] = (i == 0) ? vec[i] : sum;
	}
}


// Shift in parallel
__global__ void shiftArrayElementsInParallel(float* vec, int vecLength, float* result, short delay) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < vecLength-delay) {
        result[i] = vec[i+delay];
	}
}


// Sparse coefficients
__global__ void multiplyPolynomialCoefficientsKernel(const float* vecA, const int vecASize, const float* vecB, const int vecBSize, float* vecC) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < vecASize) {
        for (int j = 0; j < vecBSize; j++) {
			vecC[i+j] = vecC[i+j] + (vecA[i] * vecB[j]);
		}
	}
}

/**
 *	Kernel que inicia os valores de um array com 0.0.
 */
__global__ void initializeArray(float* array, int arrayLength) {
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	if (tid < arrayLength) {
		array[tid] = 0.0;
	}
}

__global__ void sumColumns(float* dev_aux1, int aux1Length, float* dev_aux_sum, int resultSize) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < resultSize && i < aux1Length) {
		dev_aux_sum[i] += dev_aux1[i];
	}
}


__global__ void sumPolynomialCoefficientsKernel(float* dev_aux1, float* dev_aux2, float* dev_aux_sum, int resultSize) {
	int i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i < resultSize) {
		dev_aux_sum[i] = dev_aux1[i] + dev_aux2[i];
	}
}

/**
 *	Não é usada.
 */
__global__ void conv(float* u, int uLength, float* v, int vLength, float* w, int wLength, int maxLength, int minLength) {
	// seguindo a referência: http://www.mathworks.com/help/techdoc/ref/conv.html
	int tid = blockDim.x * blockIdx.x + threadIdx.x;
	
	if (tid < wLength) {
		/*for (int j = 0; j < uLength; j++ ) {
			if ((tid - j) >= 0) {
				w[tid] += (u[j] * v[tid-j]);    // convolve: multiply and accumulate
			}
		}*/
		for (int j = 0; j < maxLength; j++ ) {
			if ((tid - j) >= 0 && (tid - j) < minLength) {
				if (maxLength == uLength) {
					w[tid] += (u[j] * v[tid-j]);    // convolve: multiply and accumulate
				} else {
					w[tid] += (v[j] * u[tid-j]);    // convolve: multiply and accumulate
				}				
			}
		}
	}
}


float* convSimple(float* signal, int signalLength, float* filter, int filterLength) {
    int new_size = signalLength + filterLength - 1;
	float *dev_result, *dev_signal, *dev_filter, *h_result;
	
	h_result = (float*) calloc(new_size, sizeof(float));

	cudaMalloc((void**)&dev_signal, signalLength * sizeof(float));
	cudaMalloc((void**)&dev_filter, filterLength * sizeof(float));
    cudaMalloc((void**)&dev_result, new_size * sizeof(float));
	
    // Copy host memory to device
    cudaMemcpy(dev_signal, signal, (signalLength * sizeof(float)), cudaMemcpyHostToDevice);
	cudaMemcpy(dev_filter, filter, (filterLength * sizeof(float)), cudaMemcpyHostToDevice);
	
	initializeArray<<<32,32>>>(dev_result, new_size);

	int wLength = signalLength + filterLength - 1;	
	int maxLength = (signalLength >= filterLength ? signalLength : filterLength);
	int minLength = (signalLength <= filterLength ? signalLength : filterLength);

	conv<<<32, 32>>>(dev_signal, signalLength, dev_filter, filterLength, dev_result, wLength, maxLength, minLength);
    cudaMemcpy(h_result, dev_result, (new_size * sizeof(float)), cudaMemcpyDeviceToHost);

    // cleanup memory
	cudaFree(dev_result);
    cudaFree(dev_signal);
	cudaFree(dev_filter);
	
	return h_result;
}



/**
 *	CUFFT
 **/
////////////////////////////////////////////////////////////////////////////////
//! Run a simple test for CUDA
////////////////////////////////////////////////////////////////////////////////
float* convFFT(float* signal, int signalLength, float* filter, int filterLength) {
    // Allocate host memory for the signal
    Complex* h_signal = (Complex*) malloc(sizeof(Complex) * signalLength);
	
    // Initalize the memory for the signal
    for (unsigned int i = 0; i < signalLength; ++i) {
        h_signal[i].x = signal[i];
        h_signal[i].y = 0.0;
    }

    // Allocate host memory for the filter
    Complex* h_filter_kernel = (Complex*) malloc(sizeof(Complex) * filterLength);
	
    // Initalize the memory for the filter
    for (unsigned int i = 0; i < filterLength; ++i) {
        h_filter_kernel[i].x = filter[i];
        h_filter_kernel[i].y = 0.0;
    }

    // Pad signal and filter kernel
    Complex* h_padded_signal;
    Complex* h_padded_filter_kernel;
    int new_size = PadData(h_signal, &h_padded_signal, signalLength,
                           h_filter_kernel, &h_padded_filter_kernel, filterLength);
	
    int mem_size = sizeof(Complex) * new_size;

    // Allocate device memory for signal
    Complex* d_signal;
    cudaMalloc((void**)&d_signal, mem_size);
    // Copy host memory to device
    cudaMemcpy(d_signal, h_padded_signal, mem_size,
                              cudaMemcpyHostToDevice);

    // Allocate device memory for filter kernel
    Complex* d_filter_kernel;
    cudaMalloc((void**)&d_filter_kernel, mem_size);

    // Copy host memory to device
    cudaMemcpy(d_filter_kernel, h_padded_filter_kernel, mem_size,
                              cudaMemcpyHostToDevice);

    // CUFFT plan
    cufftHandle plan;
    cufftPlan1d(&plan, new_size, CUFFT_C2C, 1);

    // Transform signal and kernel
    cufftExecC2C(plan, (cufftComplex *)d_signal, (cufftComplex *)d_signal, CUFFT_FORWARD);
    cufftExecC2C(plan, (cufftComplex *)d_filter_kernel, (cufftComplex *)d_filter_kernel, CUFFT_FORWARD);

    // Multiply the coefficients together and normalize the result
    ComplexPointwiseMulAndScale<<<32, 256>>>(d_signal, d_filter_kernel, new_size, 1.0f / new_size);

    // Check if kernel execution generated and error
    //cutilCheckMsg("Kernel execution failed [ ComplexPointwiseMulAndScale ]");

    // Transform signal back
    cufftExecC2C(plan, (cufftComplex *)d_signal, (cufftComplex *)d_signal, CUFFT_INVERSE);

    // Copy device memory to host
    Complex* h_convolved_signal = h_padded_signal;
    cudaMemcpy(h_convolved_signal, d_signal, mem_size,
                              cudaMemcpyDeviceToHost);

    // Allocate host memory for the convolution result
    Complex* h_convolved_signal_ref = (Complex*)malloc(sizeof(Complex) * signalLength);

    // Convolve on the host
	
	float* convolvedSignal = (float*) calloc(new_size, sizeof(float));
	for (int i = 0; i < new_size; i++) {
		convolvedSignal[i] = h_convolved_signal[i].x;
	}

    // check result


    //Destroy CUFFT context
    cufftDestroy(plan);

    // cleanup memory
    free(h_signal);
    free(h_filter_kernel);
    free(h_padded_signal);
    free(h_padded_filter_kernel);
    //free(h_convolved_signal_ref);
	//free(h_convolved_signal);
    cudaFree(d_signal);
    cudaFree(d_filter_kernel);
	
	return convolvedSignal;

}

// Pad data
int PadData(const Complex* signal, Complex** padded_signal, int signal_size,
            const Complex* filter_kernel, Complex** padded_filter_kernel, int filter_kernel_size) {
    int minRadius = filter_kernel_size / 2;
    int maxRadius = filter_kernel_size - minRadius;
    //int new_size = signal_size + maxRadius;
	int new_size = signal_size + filter_kernel_size - 1;
    
    // Pad signal
    Complex* new_data = (Complex*)malloc(sizeof(Complex) * new_size);
    memcpy(new_data + 0, signal, signal_size * sizeof(Complex));
    memset(new_data + signal_size, 0, (new_size - signal_size) * sizeof(Complex));
    *padded_signal = new_data;
    
    // Pad filter
    new_data = (Complex*)malloc(sizeof(Complex) * new_size);  
    /*memcpy(new_data + 0, filter_kernel + minRadius, maxRadius * sizeof(Complex));
    memset(new_data + maxRadius, 0, (new_size - filter_kernel_size) * sizeof(Complex));   
    memcpy(new_data + new_size - minRadius, filter_kernel, minRadius * sizeof(Complex));*/
	
	memcpy(new_data + 0, filter_kernel, filter_kernel_size * sizeof(Complex));
    memset(new_data + filter_kernel_size, 0, (new_size - filter_kernel_size) * sizeof(Complex));
	
    *padded_filter_kernel = new_data;
    
    return new_size;
}

////////////////////////////////////////////////////////////////////////////////
// Filtering operations
////////////////////////////////////////////////////////////////////////////////

// Computes convolution on the host - NÃO É USADA
void Convolve(const Complex* signal, int signal_size,
              const Complex* filter_kernel, int filter_kernel_size,
              Complex* filtered_signal) {
    int minRadius = filter_kernel_size / 2;
    int maxRadius = filter_kernel_size - minRadius;
    // Loop over output element indices
    for (int i = 0; i < signal_size; ++i) {
        filtered_signal[i].x = filtered_signal[i].y = 0;
        // Loop over convolution indices
        for (int j = - maxRadius + 1; j <= minRadius; ++j) {
            int k = i + j;
            if (k >= 0 && k < signal_size) 
                filtered_signal[i] = ComplexAdd(filtered_signal[i], ComplexMul(signal[k], filter_kernel[minRadius - j]));
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
// Complex operations
////////////////////////////////////////////////////////////////////////////////

// Complex addition
static __device__ __host__ inline Complex ComplexAdd(Complex a, Complex b) {
    Complex c;
    c.x = a.x + b.x;
    c.y = a.y + b.y;
    return c;
}

// Complex scale
static __device__ __host__ inline Complex ComplexScale(Complex a, float s) {
    Complex c;
    c.x = s * a.x;
    c.y = s * a.y;
    return c;
}

// Complex multiplication
static __device__ __host__ inline Complex ComplexMul(Complex a, Complex b) {
    Complex c;
    c.x = a.x * b.x - a.y * b.y;
    c.y = a.x * b.y + a.y * b.x;
    return c;
}

// Complex pointwise multiplication
static __global__ void ComplexPointwiseMulAndScale(Complex* a, const Complex* b, int size, float scale) {
    const int numThreads = blockDim.x * gridDim.x;
    const int threadID = blockIdx.x * blockDim.x + threadIdx.x;
    for (int i = threadID; i < size; i += numThreads)
        a[i] = ComplexScale(ComplexMul(a[i], b[i]), scale);     
} 



// ################## Host code

// Find Delay
/**
 Retorna o atraso (ITD) para cada ouvido usando algoritmo de busca diretamente 
 sobre a HRIR da direcao.
 D[0] = left;
 D[1] = right;
*/
short* findDelay(float** hrtf, int vecLength) {
	float** squaredElements = getSquaresOfArrayElements(hrtf[0], hrtf[1], vecLength);;

	double* et = (double*) malloc(2 * sizeof(double));
	et[0] = sumArrayElements(squaredElements[0], vecLength);
	et[1] = sumArrayElements(squaredElements[1], vecLength);
	
	float** ek = sumColumns(squaredElements[0], squaredElements[1], vecLength);;
	
	short* d = (short*) malloc(2*sizeof(short));
	for (int ear = 0; ear < 2; ear++) {	// ear
		for (int j = 0; j < vecLength; j++) {
			ek[ear][j] = ek[ear][j]/et[ear];
		}
		short* indexes = findIndexesGreaterThan(ek[ear], vecLength, 0.002f);
		d[ear] = indexes[0]-1;	
	}
	
	cleanHostMemory(ek);
	cleanHostMemory(et);
	cleanHostMemory(squaredElements);
	
	return d;
}

/**
 *	Soma paralelamente as colunas do vetor @vec de tamanho @vecLength.
 */
float** sumColumns(float* vec0, float* vec1, int vecLength) {
	// Variables
	float** h_result;

	int size = vecLength*sizeof(float);
	h_result = (float**)malloc(2 * sizeof(float*));
	h_result[0] = (float*)malloc(size);
	h_result[1] = (float*)malloc(size);
	
	for (int i = 0; i < vecLength; i++) {
		float sum0 = 0.0;
		float sum1 = 0.0;
		for (int j = 0; j < i; j++) {
			sum0 += vec0[j];
			sum1 += vec1[j];
		}		
		h_result[0][i] = (i == 0) ? vec0[i] : sum0;
		h_result[1][i] = (i == 0) ? vec1[i] : sum1;
	}
	

	/*
	float* d_vec;
	float* d_result;
	
	cudaEvent_t start, stop;
	
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	cudaEventRecord(start, 0);

	cudaMalloc((void**)&d_vec, size);
    cudaMalloc((void**)&d_result, size);

	// Copy vector from host memory to device memory
    cudaMemcpy(d_vec, vec, size, cudaMemcpyHostToDevice);

	// Invoke kernel
    sumColumnsInParallel<<<blocksPerGrid, threadsPerBlock>>>(d_vec, d_result, vecLength);

	// Copy result from device memory to host memory
    // h_result contains the result in host memory
    cudaMemcpy(h_result, d_result, size, cudaMemcpyDeviceToHost);
	
	cudaEventRecord(stop, 0);
	//cudaEventSynchronize(stop);
	
	float elapsedTime;
	cudaEventElapsedTime(&elapsedTime, start, stop);
	
	printf("\n\n -- Elapsed time: %3.15f ms --\n\n", elapsedTime);
	
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
	
	cleanDeviceMemory(d_vec);
	cleanDeviceMemory(d_result);*/
	
	return h_result;
}

/**
 *	Soma os quadrados de todos os elementos de um array.
 */
/*double sumSquaresOfArrayElements(float* vec, int vecLength) {
	float* h_result = getSquaresOfArrayElements(vec, vecLength);
	return sumArrayElements(h_result, vecLength);
}*/

/**
 *	Soma todos os elementos de um array.
 */
double sumArrayElements(float* vec, int vecLength) {
	// Sum result
	int i;
	double result = 0.0;
    for (i = 0; i < vecLength; ++i) {
		result += vec[i];
    }	
	
	return result;
}

/**
 *	Recupera os quadrados de cada elemento do array.
 */
float** getSquaresOfArrayElements(float* vec0, float* vec1, int vecLength) {
	// Variables
	float** h_result;
	int size = vecLength*sizeof(float);
	
	h_result = (float**)malloc(2 * sizeof(float*));
	h_result[0] = (float*) malloc(size);
	h_result[1] = (float*) malloc(size);
	
	for (int i = 0; i < vecLength; i++) {
		h_result[0][i] = (vec0[i] * vec0[i]);
		h_result[1][i] = (vec1[i] * vec1[i]);
	}
	

	/*
	float* d_vec0;
	float* d_vec1;
	float* d_result0;
	float* d_result1;
	
	cudaMalloc((void**)&d_vec0, size);
	cudaMalloc((void**)&d_vec1, size);
    cudaMalloc((void**)&d_result0, size);
	cudaMalloc((void**)&d_result1, size);

	// Copy vector from host memory to device memory
    cudaMemcpy(d_vec0, vec0, size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_vec1, vec1, size, cudaMemcpyHostToDevice);

	// Invoke kernel
    VecMult<<<blocksPerGrid, threadsPerBlock>>>(d_vec0, d_vec1, d_result0, d_result1, vecLength);

	// Copy result from device memory to host memory
    // h_result contains the result in host memory
    cudaMemcpy(h_result[0], d_result0, size, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_result[1], d_result1, size, cudaMemcpyDeviceToHost);
	
	cleanDeviceMemory(d_vec0);
	cleanDeviceMemory(d_vec1);
	cleanDeviceMemory(d_result0);
	cleanDeviceMemory(d_result1);*/
	
	return h_result;
}


/**
 *	Encontra os índices de um array em que os valores são maiores que @minValue.
 */
short* findIndexesGreaterThan(float* vec, int vecLength, double minValue) {
	short* indexes = (short*)malloc(vecLength*sizeof(short));
	int j = 0;
	for (short i = 0; i < vecLength; i++) {
		if (vec[i] >= minValue) {
			indexes[j++] = i;
		}
	}
	return indexes;
}


// Shift in parallel
/**
 *	Desloca os elementos de um array em paralelo.
 */
float* shiftInParallel(float* vec, int vecLength, short delay, int maxLength) {
	// Variables
	float* h_result = NULL;
	
	h_result = (float*) calloc(maxLength, sizeof(float));
	
	for (int i = 0; i < (vecLength-delay); i++) {
		h_result[i] = vec[i+delay];
	}
	
	
	/*
	float* d_vec = NULL;
	float* d_result = NULL;
	
	int vecSize = vecLength * sizeof(float);
	int size = maxLength * sizeof(float);
	
	cudaMalloc( (void**)&d_vec, vecSize );
    cudaMalloc((void**)&d_result, size);

	// Copy vector from host memory to device memory
    cudaMemcpy(d_vec, vec, vecLength*sizeof(float), cudaMemcpyHostToDevice);

	// Invoke kernel
    shiftArrayElementsInParallel<<<blocksPerGrid, threadsPerBlock>>>(d_vec, vecLength, d_result, delay);

	// Copy result from device memory to host memory
    // h_result contains the result in host memory
    cudaMemcpy(h_result, d_result, size, cudaMemcpyDeviceToHost);
	
	cleanDeviceMemory(d_vec);
	cleanDeviceMemory(d_result);*/
	
	return h_result;
}



// Sparse coefficients

/**
 *	Função que lê os coeficientes do filtro daubechies 8
 */
double** leFiltros(char* filtro, int* filtroLength) {
	double** h = NULL;
	
	if (strcmp(DB8_CONST, filtro) == 0) {
		double db8[2][8] = {
			{-0.0105, 0.0328, 0.0308, -0.1870, -0.0279, 0.6308, 0.7148, 0.2303},
			{-0.2303, 0.7148, -0.6308, -0.0279, 0.1870, 0.0308, -0.0328, -0.0105}};
	
		int channels = 2;	// número de canais nesse filtro
		int filterLength = 8;	// tamanho do maior filtro
		*filtroLength = filterLength;
		
		h = (double**) calloc(channels, sizeof(double*));
		h[0] = (double*) calloc(filterLength, sizeof(double));
		h[1] = (double*) calloc(filterLength, sizeof(double));
		
		for (int row = 0; row < channels; row++) {
			for (int col = 0; col < filterLength; col++) {
				h[row][col] = db8[row][col];
			}
		}
	}
	return h;
}

/**
 *	Função que retorna o maior elemento de um vetor
 */
int max(int* vec, int vecLength) {
	int max = 0;
	for (int k = 0; k < vecLength; k++) { 
		if (vec[k] > max) {
			max = vec[k];
		}
	}
	return max;
}

/**
 *	Função que retorna o vetor com elementos esparsos (separados por 0),
 *	de acordo com um coeficiente de esparsidade.
 */
float* spars(double* vec, int vecLength, int sparsity, int* resultLength) {
	*resultLength = (((vecLength - 1) * sparsity) + 1);
	float* y = (float*) calloc((((vecLength - 1) * sparsity) + 1) , sizeof(float));
	
	for (int i = 0; i < vecLength; i++) {
		y[(i * sparsity)] = vec[i];
	}
	
	return y;
}

/**
 *	Operação de cascateamento dos filtros utilizando a convolução
 *	com FFT implementada em CUDA.
 */
void cascataFFT(char* filtros[], int numFiltros, float** filterBank, int* filterBankLength) {
	float *convResult;

	int J = numFiltros;
	int filtroLength;
	double** h;
	int hLength;
	
	h = leFiltros(filtros[0], &hLength);
		
	float** filterBankAux = (float**) malloc((numFiltros + 1) * sizeof(float*));
	int* filterBankLengthAux = (int*) malloc((numFiltros + 1) * sizeof(int));
	filterBankAux[0] = (float*) malloc(hLength * sizeof(float));	// passa-altas
	filterBankLengthAux[0] = hLength;
	float* Gl = (float*) malloc(hLength * sizeof(float));	// passa-baixas
	
	for (int i = 0; i < hLength; i++) {
		filterBankAux[0][i] = h[1][i];
		Gl[i] = h[0][i];
	}

	float* Gl_old = NULL;
	filtroLength = hLength;
	
	cudaEvent_t	start, stop;
	float	elapsedTime;
	// start the timers
	cudaEventCreate( &start ); 
	cudaEventCreate( &stop ); 
	cudaEventRecord( start, 0 );
	
	for (int i = 1; i < numFiltros; i++) {
		Gl_old = (float*) calloc(filtroLength, sizeof(float));
		for (int j = 0; j < filtroLength; j++) {
			Gl_old[j] = Gl[j];
		}
		
		int sparsArray0Length, sparsArray1Length;
		float* sparsArray0 = spars(h[0], hLength, pow(2.0, i), &sparsArray0Length);
		float* sparsArray1 = spars(h[1], hLength, pow(2.0, i), &sparsArray1Length);
		
		int resultLength = (filtroLength + sparsArray0Length - 1);
		int resultSize = resultLength * sizeof(float);
	
		convResult = (float*) calloc(resultLength, sizeof(float));
		//convResult = convFFT(Gl_old, filtroLength, sparsArray0, sparsArray0Length);
		convResult = convSimple(Gl_old, filtroLength, sparsArray0, sparsArray0Length);
		
		filterBankAux[i] = (float*) calloc(resultLength, sizeof(float));
		//filterBankAux[i] = convFFT(Gl_old, filtroLength, sparsArray1, sparsArray1Length);
		filterBankAux[i] = convSimple(Gl_old, filtroLength, sparsArray1, sparsArray1Length);
		
		filterBankLengthAux[i] = resultLength;

		if ((i+1) != numFiltros) {
			free(Gl_old);	// Gl_old é usado fora do for após a última iteração
			free(Gl);
			
			filtroLength = resultLength;
			Gl = (float*) malloc(resultSize);
			for (int i = 0; i < resultLength; i++) {
				Gl[i] = convResult[i];
			}
			
			free(convResult);
		}
	}
		
	int sparsArray0Length;
	float* sparsArray0 = spars(h[0], hLength, pow(2.0, J-1), &sparsArray0Length);
	
	int resultLength = (filtroLength + sparsArray0Length - 1);
	
	filterBankAux[J] = (float*) calloc(resultLength, sizeof(float));
	//filterBankAux[J] = convFFT(Gl_old, filtroLength, sparsArray0, sparsArray0Length);
	filterBankAux[J] = convSimple(Gl_old, filtroLength, sparsArray0, sparsArray0Length);
	
	filterBankLengthAux[J] = resultLength;

	free(Gl_old);
	free(Gl);
	
	cudaEventRecord( stop, 0 );
	cudaEventSynchronize( stop ); 
	cudaEventElapsedTime( &elapsedTime, start, stop ); 
	//printf( "Time taken cascataFFT: %3.1f ms\n", elapsedTime );
	
	int maxLength = max(filterBankLengthAux, (numFiltros + 1));
	
	for (int i = numFiltros; i >= 0; i--) {
		filterBank[i] = (float*) calloc(maxLength, sizeof(float));
		filterBankLength[i] = filterBankLengthAux[numFiltros-i];
		for (int j = 0; j < filterBankLengthAux[numFiltros - i]; j++) {
			filterBank[i][j] = filterBankAux[numFiltros - i][j];
		}
	}
}


/**
 * Essa função computa a multiplicação das matrizes Fp e Hp. Ip_aux = Fp*Hp.
 * Resultado: [2][2][7]
 *
 *	Fp =	a	b
 *			c	d
 *
 *	Hp =	A	B
 *			C	D
 *
 *	Resultado mais rápido que método anterior com streams.
 */
float*** computeIpAux(float*** Hp, int hlength, float*** Fp, int flength) {
	float *dev_a0, *dev_b0, *dev_c0, *dev_d0, *dev_A0, *dev_B0, *dev_C0, *dev_D0; //GPU buffers for stream0 
	float *dev_aux1, *dev_aux2, *dev_aux_sum;
	
	cudaEvent_t	start, stop;
	float	elapsedTime;
	// start the timers
	cudaEventCreate( &start ); 
	cudaEventCreate( &stop ); 
	cudaEventRecord( start, 0 );
	
	int size = hlength * sizeof(float);
	int resultLength = (hlength + flength - 1);
	int resultSize = resultLength * sizeof(float);
	
	float*** Ip_aux = (float***) malloc(2 * sizeof(float**));
	Ip_aux[0] = (float**) malloc(2 * sizeof(float*));
	Ip_aux[1] = (float**) malloc(2 * sizeof(float*));
	
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			Ip_aux[i][j] = (float*) malloc(resultLength*sizeof(float*));
		}
	}
	
	cudaMalloc( (void**)&dev_a0, size );
	cudaMalloc( (void**)&dev_b0, size );
	cudaMalloc( (void**)&dev_c0, size );
	cudaMalloc( (void**)&dev_d0, size );
	cudaMalloc( (void**)&dev_A0, size );
	cudaMalloc( (void**)&dev_B0, size );
	cudaMalloc( (void**)&dev_C0, size );
	cudaMalloc( (void**)&dev_D0, size );
	
	cudaMalloc( (void**)&dev_aux1, resultSize );
	cudaMalloc( (void**)&dev_aux2, resultSize );
	cudaMalloc( (void**)&dev_aux_sum, resultSize );
	
	/*
	Ip_aux = Fp*Hp;
	O resultado é uma matriz 2x2x7.
	
	aA+bC	aB+bD
	cA+dC	cB+dD
	*/
	
/* 00 */
	cudaMemcpy(dev_a0, Fp[0][0], size, cudaMemcpyHostToDevice );
	cudaMemcpy(dev_b0, Fp[0][1], size, cudaMemcpyHostToDevice );
	cudaMemcpy(dev_c0, Fp[1][0], size, cudaMemcpyHostToDevice );
	cudaMemcpy(dev_d0, Fp[1][1], size, cudaMemcpyHostToDevice );
	cudaMemcpy(dev_A0, Hp[0][0], size, cudaMemcpyHostToDevice );
	cudaMemcpy(dev_B0, Hp[0][1], size, cudaMemcpyHostToDevice );
	cudaMemcpy(dev_C0, Hp[1][0], size, cudaMemcpyHostToDevice );
	cudaMemcpy(dev_D0, Hp[1][1], size, cudaMemcpyHostToDevice );
	
	// clean auxiliary variables
	initializeArray<<<4,4>>>(dev_aux1, resultLength);
	initializeArray<<<4,4>>>(dev_aux2, resultLength);
	
	// Invoke kernel
    multiplyPolynomialCoefficientsKernel<<<hlength, 1>>>(dev_a0, hlength, dev_A0, hlength, dev_aux1);
    multiplyPolynomialCoefficientsKernel<<<hlength, 1>>>(dev_b0, hlength, dev_C0, hlength, dev_aux2);	
	
	sumPolynomialCoefficientsKernel<<<resultSize, 1>>>(dev_aux1, dev_aux2, dev_aux_sum, resultLength);
	
	// copy the data from device to locked memory
	cudaMemcpy(Ip_aux[0][0], dev_aux_sum, resultSize, cudaMemcpyDeviceToHost);	


/* 01 */
	// clean auxiliary variables
	initializeArray<<<4,4>>>(dev_aux1, resultLength);
	initializeArray<<<4,4>>>(dev_aux2, resultLength);

	// Invoke kernel
    multiplyPolynomialCoefficientsKernel<<<hlength, 1>>>(dev_a0, hlength, dev_B0, hlength, dev_aux1);
    multiplyPolynomialCoefficientsKernel<<<hlength, 1>>>(dev_b0, hlength, dev_D0, hlength, dev_aux2);
	
	sumPolynomialCoefficientsKernel<<<resultSize, 1>>>(dev_aux1, dev_aux2, dev_aux_sum, resultLength);
	
	// copy the data from device to locked memory
	cudaMemcpy(Ip_aux[0][1], dev_aux_sum, resultSize, cudaMemcpyDeviceToHost);


/* 10 */
	// clean auxiliary variables
	initializeArray<<<4,4>>>(dev_aux1, resultLength);
	initializeArray<<<4,4>>>(dev_aux2, resultLength);
	
	// Invoke kernel
    multiplyPolynomialCoefficientsKernel<<<hlength, 1>>>(dev_c0, hlength, dev_A0, hlength, dev_aux1);
    multiplyPolynomialCoefficientsKernel<<<hlength, 1>>>(dev_d0, hlength, dev_C0, hlength, dev_aux2);	
	
	sumPolynomialCoefficientsKernel<<<resultSize, 1>>>(dev_aux1, dev_aux2, dev_aux_sum, resultLength);
	
	// copy the data from device to locked memory
	cudaMemcpy(Ip_aux[1][0], dev_aux_sum, resultSize, cudaMemcpyDeviceToHost);	


/* 11 */
	// clean auxiliary variables
	initializeArray<<<4,4>>>(dev_aux1, resultLength);
	initializeArray<<<4,4>>>(dev_aux2, resultLength);

	// Invoke kernel
    multiplyPolynomialCoefficientsKernel<<<hlength, 1>>>(dev_c0, hlength, dev_B0, hlength, dev_aux1);
    multiplyPolynomialCoefficientsKernel<<<hlength, 1>>>(dev_d0, hlength, dev_D0, hlength, dev_aux2);	
	
	sumPolynomialCoefficientsKernel<<<resultSize, 1>>>(dev_aux1, dev_aux2, dev_aux_sum, resultLength);
	
	// copy the data from device to locked memory
	cudaMemcpy(Ip_aux[1][1], dev_aux_sum, resultSize, cudaMemcpyDeviceToHost);	

	cudaEventRecord( stop, 0 );
	cudaEventSynchronize( stop ); 
	cudaEventElapsedTime( &elapsedTime, start, stop ); 
	//printf( "Time taken computeIpAux: %3.1f ms\n", elapsedTime );
	
	cudaFree(dev_a0);
	cudaFree(dev_b0);
	cudaFree(dev_c0);
	cudaFree(dev_d0);
	cudaFree(dev_A0);
	cudaFree(dev_B0);
	cudaFree(dev_C0);
	cudaFree(dev_D0);
	cudaFree(dev_aux1);
	cudaFree(dev_aux2);
	cudaFree(dev_aux_sum);	
	
	//printf("Error computeIpAux: %s\n\n", cudaGetErrorString(cudaGetLastError()));
	
	return Ip_aux;
}


/**
 * Essa função computa a multiplicação das matrizes Rp, Fp e Ip_aux. Ip_aux = Fp*Hp.
 * Resultado: [2][237]
 *
 *	Rp =	X
 *			Y
 *
 *	Fp =	a	b
 *			c	d
 *
 *	Ip =	A	B
 *			C	D
 *
 */
float** calculateGp(float** Rp, int rpLength, float*** Fp, int fpLength, float*** Ip, int ipLength) {
	// Rp
	float *dev_x0 = NULL;
	float *dev_y0 = NULL;

    // Fp
	float *dev_a0 = NULL;
	float *dev_b0 = NULL;
	float *dev_c0 = NULL;
	float *dev_d0 = NULL;
	
	// Ip
	float *dev_A0 = NULL;
	float *dev_B0 = NULL;
	float *dev_C0 = NULL;
	float *dev_D0 = NULL;
	
	float *dev_aux1, *dev_aux2, *dev_aux3, *dev_aux4, *dev_aux_sum1, *dev_aux_sum2 ;
	float **finalResult;
	
	cudaEvent_t	start, stop;
	float	elapsedTime;
	// start the timers
	cudaEventCreate( &start ); 
	cudaEventCreate( &stop ); 
	cudaEventRecord( start, 0 );

	int lengthX = (rpLength % 2 == 0 ? rpLength/2 : (rpLength/2+1));
	int lengthY = rpLength/2;
	int sizeX = lengthX * sizeof(float);
	int sizeY = lengthY * sizeof(float);
	int sizeFpElement = fpLength * sizeof(float);
	int sizeIpElement = ipLength * sizeof(float);
	int partialResultLengthX = lengthX + fpLength - 1;
	int partialResultSizeX = partialResultLengthX * sizeof(float);
	int partialResultLengthY = lengthY + fpLength - 1;
	int partialResultSizeY = partialResultLengthY * sizeof(float);
	int finalResultLength = (partialResultLengthX + ipLength -1);
	int finalResultSize = finalResultLength * sizeof(float);
	
	cudaMalloc( (void**)&dev_x0, sizeX );
	cudaMalloc( (void**)&dev_y0, sizeY );
	cudaMalloc( (void**)&dev_a0, sizeFpElement );
	cudaMalloc( (void**)&dev_b0, sizeFpElement );
	cudaMalloc( (void**)&dev_c0, sizeFpElement );
	cudaMalloc( (void**)&dev_d0, sizeFpElement );
	cudaMalloc( (void**)&dev_A0, sizeIpElement );
	cudaMalloc( (void**)&dev_B0, sizeIpElement );
	cudaMalloc( (void**)&dev_C0, sizeIpElement );
	cudaMalloc( (void**)&dev_D0, sizeIpElement );
	
	cudaMalloc( (void**)&dev_aux1, partialResultSizeX );
	cudaMalloc( (void**)&dev_aux2, partialResultSizeY );
	cudaMalloc( (void**)&dev_aux3, partialResultSizeX );
	cudaMalloc( (void**)&dev_aux4, partialResultSizeY );
	
	cudaMalloc( (void**)&dev_aux_sum1, partialResultSizeX );
	cudaMalloc( (void**)&dev_aux_sum2, partialResultSizeX );
	
	// Inicialização de array no dispositivo para evitar erros em leituras
	initializeArray<<<16,16>>>(dev_aux1, partialResultLengthX);
	initializeArray<<<16,16>>>(dev_aux2, partialResultLengthX);
	initializeArray<<<16,16>>>(dev_aux3, partialResultLengthX);
	initializeArray<<<16,16>>>(dev_aux4, partialResultLengthX);
	initializeArray<<<16,16>>>(dev_aux_sum1, partialResultLengthX);
	initializeArray<<<16,16>>>(dev_aux_sum2, partialResultLengthX);
	
	/*
	Gp = Rp*Fp*Ip;
	O resultado é uma matriz 2x237.
	
	----- Parte 1
	  Rp	  Fp
	(x	y)	(a	c)
			(b	d)
	
	Rp*Fp
		xa+yb	xc+yd
	
	*/
	cudaMemcpy(dev_x0, Rp[0], sizeX, cudaMemcpyHostToDevice );
	cudaMemcpy(dev_y0, Rp[1], sizeY, cudaMemcpyHostToDevice );
	
	cudaMemcpy(dev_a0, Fp[0][0], sizeFpElement, cudaMemcpyHostToDevice );
	cudaMemcpy(dev_b0, Fp[0][1], sizeFpElement, cudaMemcpyHostToDevice );
	cudaMemcpy(dev_c0, Fp[1][0], sizeFpElement, cudaMemcpyHostToDevice );
	cudaMemcpy(dev_d0, Fp[1][1], sizeFpElement, cudaMemcpyHostToDevice );
	
	cudaMemcpy(dev_A0, Ip[0][0], sizeIpElement, cudaMemcpyHostToDevice );
	cudaMemcpy(dev_B0, Ip[0][1], sizeIpElement, cudaMemcpyHostToDevice );
	cudaMemcpy(dev_C0, Ip[1][0], sizeIpElement, cudaMemcpyHostToDevice );
	cudaMemcpy(dev_D0, Ip[1][1], sizeIpElement, cudaMemcpyHostToDevice );
	
	// Invoke kernel
    multiplyPolynomialCoefficientsKernel<<<16, 16>>>(dev_x0, lengthX, dev_a0, fpLength, dev_aux1);	
    multiplyPolynomialCoefficientsKernel<<<16, 16>>>(dev_y0, lengthY, dev_b0, fpLength, dev_aux2);
	multiplyPolynomialCoefficientsKernel<<<16, 16>>>(dev_x0, lengthX, dev_c0, fpLength, dev_aux3);
    multiplyPolynomialCoefficientsKernel<<<16, 16>>>(dev_y0, lengthY, dev_d0, fpLength, dev_aux4);
	
	sumPolynomialCoefficientsKernel<<<16, 16>>>(dev_aux1, dev_aux2, dev_aux_sum1, partialResultLengthX);
	sumPolynomialCoefficientsKernel<<<16, 16>>>(dev_aux3, dev_aux4, dev_aux_sum2, partialResultLengthX);
		
	// Re-Inicialização de array no dispositivo para evitar erros em leituras
	cudaFree(dev_aux1);
	cudaFree(dev_aux2);
	cudaFree(dev_aux3);
	cudaFree(dev_aux4);
	
	cudaMalloc( (void**)&dev_aux1, finalResultSize );
	cudaMalloc( (void**)&dev_aux2, finalResultSize );
	cudaMalloc( (void**)&dev_aux3, finalResultSize );
	cudaMalloc( (void**)&dev_aux4, finalResultSize );
	
	initializeArray<<<16,16>>>(dev_aux1, finalResultLength);
	initializeArray<<<16,16>>>(dev_aux2, finalResultLength);
	initializeArray<<<16,16>>>(dev_aux3, finalResultLength);
	initializeArray<<<16,16>>>(dev_aux4, finalResultLength);
	
	/*
	Gp = Rp*Fp*Ip;
	O resultado é uma matriz 2x237.
	
	----- Parte 1
	  Rp	  Fp
	(x	y)	(a	c)
			(b	d)
	
	Rp*Fp
		xa+yb	xc+yd
		  u		  v
	
	----- Parte 2	  	  
	Rp*Fp	  Hp
	(u	v)	(A	C)
			(B	D)
	
	Rp*Fp*Ip
		uA+vB	uC+vD
	
	*/
	multiplyPolynomialCoefficientsKernel<<<16, 16>>>(dev_aux_sum1, partialResultLengthX, dev_A0, ipLength, dev_aux1);	
    multiplyPolynomialCoefficientsKernel<<<16, 16>>>(dev_aux_sum2, partialResultLengthX, dev_B0, ipLength, dev_aux2);
	multiplyPolynomialCoefficientsKernel<<<16, 16>>>(dev_aux_sum1, partialResultLengthX, dev_C0, ipLength, dev_aux3);
    multiplyPolynomialCoefficientsKernel<<<16, 16>>>(dev_aux_sum2, partialResultLengthX, dev_D0, ipLength, dev_aux4);
	
	// Cleaning data of auxiliary vectors
	cudaFree(dev_aux_sum1);
	cudaFree(dev_aux_sum2);
	cudaMalloc( (void**)&dev_aux_sum1, finalResultSize );
	cudaMalloc( (void**)&dev_aux_sum2, finalResultSize );
	
	// Reinicializa arrays u e v (resultado de Rp*Fp)
	initializeArray<<<16,16>>>(dev_aux_sum1, finalResultLength);
	initializeArray<<<16,16>>>(dev_aux_sum2, finalResultLength);
	
	sumPolynomialCoefficientsKernel<<<16, 16>>>(dev_aux1, dev_aux2, dev_aux_sum1, finalResultLength);
	sumPolynomialCoefficientsKernel<<<16, 16>>>(dev_aux3, dev_aux4, dev_aux_sum2, finalResultLength);
	
	finalResult = (float**) malloc(2 * sizeof(float*));
	finalResult[0] = (float*) malloc(finalResultSize);
	finalResult[1] = (float*) malloc(finalResultSize);
	
	// copy the data from device to locked memory
	cudaMemcpy(finalResult[0], dev_aux_sum1, finalResultSize, cudaMemcpyDeviceToHost);
	cudaMemcpy(finalResult[1], dev_aux_sum2, finalResultSize, cudaMemcpyDeviceToHost);
	
	cudaEventRecord( stop, 0 );
	cudaEventSynchronize( stop ); 
	cudaEventElapsedTime( &elapsedTime, start, stop ); 
	//printf( "Time taken calculateGp: %3.1f ms\n", elapsedTime );
		
	cudaFree(dev_x0);
	cudaFree(dev_y0);
	cudaFree(dev_a0);
	cudaFree(dev_b0);
	cudaFree(dev_c0);
	cudaFree(dev_d0);
	cudaFree(dev_A0);
	cudaFree(dev_B0);
	cudaFree(dev_C0);
	cudaFree(dev_D0);
	
	cudaFree(dev_aux1);
	cudaFree(dev_aux2);
	cudaFree(dev_aux3);
	cudaFree(dev_aux4);
	cudaFree(dev_aux_sum1);
	cudaFree(dev_aux_sum2);

	//printf("Error calculateGp: %s\n\n", cudaGetErrorString(cudaGetLastError()));
	
	return finalResult;
}

/**
 *	Função que realiza a operação de decimação polinomial.
 *
 *  Retorna os coeficientes G (matriz 2xK) que correspondem ao 
 *	sistema R decomposto pelo banco de 2 canais H.
 */
void dec_poly(int hlength, float* sist, int sistLength, float** G, int* G_size, float*** Hp, float*** Fp) {	
	float** Rp = (float**) calloc(2, sizeof(float*));
	Rp[0] = (float*) calloc((sistLength % 2 == 0 ? (sistLength/2) : (sistLength/2 + 1)), sizeof(float));
	Rp[1] = (float*) calloc((sistLength/2), sizeof(float));

	for (int i = 0; i < sistLength; i++) {
		if (i % 2 == 0) {
			Rp[0][i/2] = sist[i];
		} else {
			Rp[1][i/2] = sist[i];
		}
	}
	
	/*
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% decomposicao polifasica do sistema R - matriz Rp
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	*/
	/*for (int i = 0; i < 10; i++) {
		printf("sist[%d] = %1.15f\n", i, sist[i]);
	}*/
	int m = 2;
	int n = hlength;

	/*
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% decomposicao polifasica dos bancos de sintese e analise
	% matrizes Fp e Hp
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	*/

	/*
	Ip_aux = Fp*Hp;
	O resultado é uma matriz 2x2x7.
	
	H[0][1]*F[0][0] + H[0][0]*F[1][0]	H[0][1]*F[0][1] + H[0][0]*F[1][1]
	H[1][1]*F[0][0] + H[1][0]*F[1][0]	H[1][1]*F[0][1] + H[1][0]*F[1][1]
	*/
	// Hp e Fp são matrizes (cubos) de dimensão: [2][2][7]	
	//float*** ip_aux_streamed = computeIpAuxStream(Hp, hlength/2, Fp, hlength/2);
	
	INIT_VARIABLES;
	INIT_RUNTIME;
	float*** ip_aux = computeIpAux(Hp, hlength/2, Fp, hlength/2);
	END_RUNTIME;
	printf("\n[computeIpAux]: ");
	PRINT_RUNTIME;
	
	/*
	for (int dim1 = 0; dim1 < 2; dim1++) {
		for (int dim2 = 0; dim2 < 2; dim2++) {
			printf("\n\n[%d][%d]---------------------\n", dim1, dim2);
			for (int i = 0; i < 7; i++) {
				printf("ip_aux_streamed[%d] = %1.15f\t\tIp_aux[%d] = %1.15f\n", i, ip_aux_streamed[dim1][dim2][i], i, ip_aux[dim1][dim2][i]);
			}
		}
	}
	printf("\n\n");
	*/
	
	// Gp = Rp*Fp*Ip
	INIT_RUNTIME;
	float** Gp = calculateGp(Rp, sistLength, Fp, (hlength/2), ip_aux, (hlength-1));
	END_RUNTIME;
	printf("\n[calculateGp]: ");
	PRINT_RUNTIME;
	
	int atraso = round((n-m)/2.0);
	double esparsidade = 2.0;
	
	*G_size = ceil( ((sistLength + atraso +1)/esparsidade) +1);
	
	G[0] = (float*) malloc(*G_size * sizeof(float));
	G[1] = (float*) malloc(*G_size * sizeof(float));
	
	for (int i = 0; i < *G_size; i++) {
		G[0][i] = Gp[0][i+atraso];
		G[1][i] = Gp[1][i+atraso];
	}
}

/**
 *	Obtém os coeficientes esparsos que equivalem o sistema ho1d.
 */
void coef_spars(char* filtro[], int numFiltros, float* ho1d, int ho1dLength, float** G_aux, int* G_size) {
	int mesmofiltro = 0;
	float* sist = ho1d;
	double** h;
	float** G;
	float** coefSpars;
	int filtroLength;
	int sistLength = ho1dLength;
	
	coefSpars = (float**) malloc((numFiltros+1) * sizeof(float*));
	
	G = (float**) malloc(2 * sizeof(float*));	// sempre tem tamanho 2
	
	if (!mesmofiltro) {
		h = leFiltros(filtro[0], &filtroLength);
		mesmofiltro = (0 < numFiltros-1 && strcmp(filtro[0], filtro[1]) == 0);
	}
	
	// Banco de análise
	double** H = (double**) malloc(2 * sizeof(double*));
	H[0] = (double*) malloc(filtroLength * sizeof(double));
	H[1] = (double*) malloc(filtroLength * sizeof(double));

	// Banco de síntese
	double** F = (double**) malloc(2 * sizeof(double*));
	F[0] = (double*) malloc(filtroLength * sizeof(double));
	F[1] = (double*) malloc(filtroLength * sizeof(double));
	
	/*
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Obtencao do banco de filtros
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	*/
	int power = 1;
	for (int i = filtroLength-1; i >= 0; i--) {
		int reverseIndex = filtroLength-i-1;
		double powerValue = pow((double)-1, (double)power++);
		if (h[0][reverseIndex] + h[1][i]*powerValue != 0) {
			break;	// Não é Daubechies
		} else {
			H[0][reverseIndex] = h[0][reverseIndex];
			H[1][reverseIndex] = h[0][i]*powerValue;
			F[0][reverseIndex] = h[1][reverseIndex]*powerValue;
			F[1][reverseIndex] = h[1][i];
		}
	}
	
	float*** Hp = (float***) malloc(2 * sizeof(float**));
	Hp[0] = (float**) malloc(2 * sizeof(float*));
	Hp[1] = (float**) malloc(2 * sizeof(float*));
	
	float*** Fp = (float***) malloc(2 * sizeof(float**));
	Fp[0] = (float**) malloc(2 * sizeof(float*));
	Fp[1] = (float**) malloc(2 * sizeof(float*));
	
	for (int i = 0; i < 2; i++) {
		Hp[0][i] = (float*) malloc((filtroLength/2) * sizeof(float));
		Hp[1][i] = (float*) malloc((filtroLength/2) * sizeof(float));
		Fp[0][i] = (float*) malloc((filtroLength/2) * sizeof(float));
		Fp[1][i] = (float*) malloc((filtroLength/2) * sizeof(float));
		for (int j = 0; j < (filtroLength/2); j++) {
			Hp[0][i][j] = H[i][j*2];
			Hp[1][i][j] = H[i][j*2+1];
			
			if (i == 0) {
				Fp[0][i][j] = F[0][j*2+1];
				Fp[1][i][j] = F[1][j*2+1];
			} else {
				Fp[0][i][j] = (-1) * F[1][filtroLength*i - (j*2+1)];
				Fp[1][i][j] = F[0][filtroLength*i - (j*2+1)];
			}
		}
	}
	
	INIT_VARIABLES;
	INIT_RUNTIME;
	
	for (int j = 0; j < numFiltros; j++) {		
		dec_poly(filtroLength, sist, sistLength, G, &G_size[numFiltros-j], Hp, Fp);
		
		sistLength = G_size[numFiltros-j];
		coefSpars[numFiltros-j] = G[1];
		
		sist = NULL;
		sist = G[0];
		 
		// Sempre usando o mesmo filtro
		/*if (!mesmofiltro) {
			h = leFiltros(filtro[j], &filtroLength);
			mesmofiltro = (j < numFiltros-1 && strcmp(filtro[j], filtro[j+1]) == 0);
		}*/	
	}
	
	END_RUNTIME;
	printf("\n[dec_poly]: ");
	PRINT_RUNTIME;
	
	G_size[0] = G_size[1];
	coefSpars[0] = (float*) malloc(G_size[0] * sizeof(float));
	for (int k = 0; k < G_size[0]; k++) {
		coefSpars[0][k] = sist[k];
	}
	
	int maxG_size = max(G_size, (numFiltros+1));
	for (int k = 0; k < (numFiltros+1); k++) {
		G_aux[k] = (float*) malloc(maxG_size * sizeof(float));
		int cont = 0;
		while (cont < G_size[k]) {
			G_aux[k][cont] = coefSpars[k][cont];
			cont++;
		}
		while (cont < maxG_size) {
			G_aux[k][cont] = 0.0;
			cont++;
		}
	}	
}

/**
 *	Função que retorna os primeiros X elementos de um vetor.
 */
double* getPartialVector(double* vec, int numOfElements) {
    double* partialVec = (double*) malloc(numOfElements * sizeof(double));
	for (int i = 0; i < numOfElements; i++) {
		partialVec[i] = vec[i];
	}
	return partialVec;
}

/**
 *	Função que retorna a resposta impulsiva a partir dos coeficientes
 *	esparsos.
 */
float* resp_imp(char* filtros[], int numFiltros, double** G, int* G_size, int* resultLength) {
	float** filterBank = (float**) calloc((numFiltros + 1), sizeof(float));
	int* filterBankLength = (int*) calloc((numFiltros + 1), sizeof(int));
		
	cascataFFT(filtros, numFiltros, filterBank, filterBankLength);
	
	// TODO implementar calc_delta
	int atrasos[5] = {1, 1, 8, 22, 50};
	
	int* L = (int*) calloc(numFiltros + 1, sizeof(int));
	L[0] = (int) pow(2.0, numFiltros);
	for (int i = 1; i < (numFiltros + 1); i++) {
		L[i] = (int) pow(2.0, (numFiltros+1)-i);
	}
	
	float** gaux = (float**) calloc((numFiltros+1), sizeof(float*));
	float** r = (float**) calloc((numFiltros+1), sizeof(float*));
	int* r_sizes = (int*) calloc((numFiltros+1), sizeof(int));
	
	for (int i = 0; i < numFiltros+1; i++) {
		int partialVecLength = (G_size[i]+atrasos[i]);
		int resultLength;
		gaux[i] = spars(getPartialVector(G[i], partialVecLength ), partialVecLength, L[i], &resultLength);
		
		// conv(filterBank[i], gaux[i])	
		int convLength = (filterBankLength[i] + resultLength - 1);
	
		//r[i] = convFFT(filterBank[i], filterBankLength[i], gaux[i], resultLength);
		r[i] = convSimple(filterBank[i], filterBankLength[i], gaux[i], resultLength);
		
		r_sizes[i] = convLength;	
	}

	int maxR = max(r_sizes, numFiltros+1);
	float* res = (float*) calloc(maxR, sizeof(float));
	float** raux = (float**) malloc((numFiltros + 1) * sizeof(float*));
	
	for (int i = 0; i < (numFiltros + 1); i++) {
		raux[i] = (float*) malloc(maxR * sizeof(float));
		for (int j = 0; j < maxR; j++) {
			if (j < r_sizes[i]) {
				raux[i][j] = r[i][j];
			} else {
				raux[i][j] = 0.0;
			}
		}
	}
	
	float* dev_r;
	float* dev_aux_sum;
	
	cudaMalloc( (void**)&dev_aux_sum, maxR * sizeof(float) );
	initializeArray<<<32,ceil(maxR/32.0)>>>(dev_aux_sum, maxR);
	
	for (int i = 0; i < (numFiltros + 1); i++) {
		cudaMalloc( (void**)&dev_r, r_sizes[i] * sizeof(float) );
		initializeArray<<<32,ceil(maxR/32.0)>>>(dev_r, r_sizes[i]);
		cudaMemcpy(dev_r, r[i], r_sizes[i] * sizeof(float), cudaMemcpyHostToDevice );
		sumColumns<<<32,ceil(maxR/32.0)>>>(dev_r, r_sizes[i], dev_aux_sum, maxR);
		cudaFree(dev_r);
	}
		
	cudaMemcpy(res, dev_aux_sum, maxR * sizeof(float), cudaMemcpyDeviceToHost);

	cudaFree(dev_aux_sum);	
	cudaFree(dev_r);
	
	int maxHlength = max(filterBankLength, numFiltros+1);
	float* resFinal = (float*) calloc((maxR - maxHlength), sizeof(float));
	
	for (int i = 0; i < (maxR-maxHlength); i++) {
		resFinal[i] = res[maxHlength + i];
	}
	
	*resultLength = (maxR-maxHlength);
	
	return resFinal;
}




// Memory management
/**
 *	Libera memória do Host.
 */
void cleanHostMemory(void* h) {
    // Free host memory
    if (h) {
        free(h);
	}
}

/**
 *	Libera memória do Device.
 */
void cleanDeviceMemory(void* d) {
    // Free device memory
    if (d) {
        cudaFree(d);
	}
	
    cudaThreadExit();	
}


// CUDA Error
void checkCUDAError(const char *msg) {
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err)  {
        fprintf(stderr, "Cuda error: %s: %s.\n", msg, cudaGetErrorString( err) );
        exit(EXIT_FAILURE);
    }                         
}










/*
 *	NÃO ESTÁ SENDO USADO
 *
 */
float* multiplyPolynomialCoefficients(float* vecA, const int vecASize, float* vecB, const int vecBSize) {
	int vecCSize = vecASize + vecBSize - 1;
	int size = vecCSize * sizeof(float);
	
	float* vecC = (float*) malloc(size);
	
	for (int i = 0; i < vecASize; i++) {
		for (int j = 0; j < vecBSize; j++) {
			vecC[i+j] = vecC[i+j] + (vecA[i] * vecB[j]);
		}
	}

	return vecC;
}

/**
 *	Operação de convolução executada no HOST.
 */
float* convHost(float* Gl_old, int filtroLength, float* sparsArray1, int resultLength) {
	float* convResult = (float*) calloc(resultLength, sizeof(float));
	for (int k = 0; k < resultLength; k++ ) {
		convResult[k] = 0;                 // set to zero before sum
		for (int j = 0; j < filtroLength; j++ ) {
			if ((k - j) >= 0) {
				convResult[k] += (Gl_old[j] * sparsArray1[k-j]);    // convolve: multiply and accumulate
			}
		}
	}
	return convResult;
}

/**
 *	Operação de cascateamento dos filtros usando CUDA e a operação de
 *	convolução em paralelo:
 *  __global__ void conv(float* u, int uLength, float* v, int vLength, float* w, int wLength);
 */
/*void cascata(char* filtros[], int numFiltros, float** filterBank, int* filterBankLength) {
	float *dev_filtro, *dev_spars_array0, *dev_spars_array1, *dev_result0, *dev_result1;
	float *convResult;

	int J = numFiltros;
	int filtroLength;
	double** h;
	int hLength;
	
	h = leFiltros(filtros[0], &hLength);
		
	float** filterBankAux = (float**) malloc((numFiltros + 1) * sizeof(float*));
	int* filterBankLengthAux = (int*) malloc((numFiltros + 1) * sizeof(int));
	filterBankAux[0] = (float*) malloc(hLength * sizeof(float));	// passa-altas
	filterBankLengthAux[0] = hLength;
	float* Gl = (float*) malloc(hLength * sizeof(float));	// passa-baixas
	
	for (int i = 0; i < hLength; i++) {
		filterBankAux[0][i] = h[1][i];
		Gl[i] = h[0][i];
	}

	float* Gl_old = NULL;
	filtroLength = hLength;
	
	cudaEvent_t	start, stop;
	float	elapsedTime;
	// start the timers
	cudaEventCreate( &start ); 
	cudaEventCreate( &stop ); 
	cudaEventRecord( start, 0 );
	
	for (int i = 1; i < numFiltros; i++) {
		Gl_old = (float*) calloc(filtroLength, sizeof(float));
		for (int j = 0; j < filtroLength; j++) {
			Gl_old[j] = Gl[j];
		}
		
		int sparsArray0Length, sparsArray1Length;
		float* sparsArray0 = spars(h[0], hLength, pow(2.0, i), &sparsArray0Length);
		float* sparsArray1 = spars(h[1], hLength, pow(2.0, i), &sparsArray1Length);
		
		int filtroSize = filtroLength * sizeof(float);
		int sparsArray0Size = sparsArray0Length * sizeof(float);
		int sparsArray1Size = sparsArray1Length * sizeof(float);
		int resultLength = (filtroLength + sparsArray0Length - 1);
		int resultSize = resultLength * sizeof(float);
	
		cudaMalloc( (void**)&dev_filtro, filtroSize );
		cudaMalloc( (void**)&dev_spars_array0, sparsArray0Size );
		cudaMalloc( (void**)&dev_spars_array1, sparsArray1Size );
		cudaMalloc( (void**)&dev_result0, resultSize );
		cudaMalloc( (void**)&dev_result1, resultSize );
	
		cudaMemcpy(dev_filtro, Gl_old, filtroSize, cudaMemcpyHostToDevice );
		cudaMemcpy(dev_spars_array0, sparsArray0, sparsArray0Size, cudaMemcpyHostToDevice );
		cudaMemcpy(dev_spars_array1, sparsArray1, sparsArray1Size, cudaMemcpyHostToDevice );
		
		initializeArray<<<16,16>>>(dev_result0, resultLength);
		initializeArray<<<16,16>>>(dev_result1, resultLength);
		
		conv<<<16,16>>>(dev_filtro, filtroLength, dev_spars_array0, sparsArray0Length, dev_result0, resultLength);
		conv<<<16,16>>>(dev_filtro, filtroLength, dev_spars_array1, sparsArray1Length, dev_result1, resultLength);
		
		convResult = (float*) calloc(resultLength, sizeof(float));
		filterBankAux[i] = (float*) calloc(resultLength, sizeof(float));
		cudaMemcpy(convResult, dev_result0, resultSize, cudaMemcpyDeviceToHost);
		cudaMemcpy(filterBankAux[i], dev_result1, resultSize, cudaMemcpyDeviceToHost);
		filterBankLengthAux[i] = resultLength;

		if ((i+1) != numFiltros) {
			free(Gl_old);	// Gl_old é usado fora do for após a última iteração
			free(Gl);
			
			filtroLength = resultLength;
			Gl = (float*) malloc(resultSize);
			for (int i = 0; i < resultLength; i++) {
				Gl[i] = convResult[i];
			}
			
			free(convResult);
		}
		
		cudaFree(dev_filtro);
		cudaFree(dev_spars_array0);
		cudaFree(dev_spars_array1);
		cudaFree(dev_result0);
		cudaFree(dev_result1);
	}
		
	int sparsArray0Length;
	float* sparsArray0 = spars(h[0], hLength, pow(2.0, J-1), &sparsArray0Length);
		
	int filtroSize = filtroLength * sizeof(float);
	int sparsArray0Size = sparsArray0Length * sizeof(float);
	int resultLength = (filtroLength + sparsArray0Length - 1);
	int resultSize = resultLength * sizeof(float);
	
	cudaMalloc( (void**)&dev_filtro, filtroSize );
	cudaMalloc( (void**)&dev_spars_array0, sparsArray0Size );
	cudaMalloc( (void**)&dev_result0, resultSize );
	
	cudaMemcpy(dev_filtro, Gl_old, filtroSize, cudaMemcpyHostToDevice );
	cudaMemcpy(dev_spars_array0, sparsArray0, sparsArray0Size, cudaMemcpyHostToDevice );
		
	initializeArray<<<16,16>>>(dev_result0, resultLength);
		
	conv<<<16,16>>>(dev_filtro, filtroLength, dev_spars_array0, sparsArray0Length, dev_result0, resultLength);
	
	filterBankAux[J] = (float*) calloc(resultLength, sizeof(float));
	cudaMemcpy(filterBankAux[J], dev_result0, resultSize, cudaMemcpyDeviceToHost);
	filterBankLengthAux[J] = resultLength;

	free(Gl_old);
	free(Gl);
		
	cudaFree(dev_filtro);
	cudaFree(dev_spars_array0);
	cudaFree(dev_result0);
	
	cudaEventRecord( stop, 0 );
	cudaEventSynchronize( stop ); 
	cudaEventElapsedTime( &elapsedTime, start, stop ); 
	//printf( "Time taken cascata: %3.1f ms\n", elapsedTime );
	
	int maxLength = max(filterBankLengthAux, (numFiltros + 1));
	
	for (int i = numFiltros; i >= 0; i--) {
		filterBank[i] = (float*) calloc(maxLength, sizeof(float));
		filterBankLength[i] = filterBankLengthAux[numFiltros-i];
		for (int j = 0; j < filterBankLengthAux[numFiltros - i]; j++) {
			filterBank[i][j] = filterBankAux[numFiltros - i][j];
		}
	}
}/*

/**
 * Essa função computa a multiplicação das matrizes Fp e Hp. Ip_aux = Fp*Hp.
 * Resultado: [2][2][7]
 *
 *	Fp =	a	b
 *			c	d
 *
 *	Hp =	A	B
 *			C	D
 *
 *
 *	NÃO ESTÁ SENDO USADO
 */
float*** computeIpAuxStream(float*** Hp, int hlength, float*** Fp, int flength) {
	float *host_A = NULL;
	float *host_B = NULL;
	float *host_C = NULL;
	float *host_D = NULL;
	float *host_a = NULL;
	float *host_b = NULL;
	float *host_c = NULL;
	float *host_d = NULL;
	float *dev_a0, *dev_b0, *dev_c0, *dev_d0, *dev_A0, *dev_B0, *dev_C0, *dev_D0; //GPU buffers for stream0 
	float *dev_aux1, *dev_aux2, *dev_aux_sum1, *dev_aux_sum2, *dev_aux_sum3, *dev_aux_sum4;
	
	cudaDeviceProp prop; 
	int whichDevice; 
	cudaGetDevice( &whichDevice ); 
	cudaGetDeviceProperties( &prop, whichDevice );
	
	if (!prop.deviceOverlap) { 
		//printf( "Device will not handle overlaps, so no speed up from streams\n" );
		return NULL;
	}
	
	cudaEvent_t	start, stop;
	float	elapsedTime;
	// start the timers
	cudaEventCreate( &start ); 
	cudaEventCreate( &stop ); 
	cudaEventRecord( start, 0 );
	
	// initialize the streams
	cudaStream_t stream0, stream1; 
	cudaStreamCreate( &stream0 );
	cudaStreamCreate( &stream1 );
	
	int size = hlength * sizeof(float);
	int resultLength = (hlength + flength - 1);
	int resultSize = resultLength * sizeof(float);
	
	float*** Ip_aux = (float***) malloc(2 * sizeof(float**));
	Ip_aux[0] = (float**) malloc(2 * sizeof(float*));
	Ip_aux[1] = (float**) malloc(2 * sizeof(float*));
	
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			Ip_aux[i][j] = (float*) malloc(resultLength*sizeof(float*));
		}
	}
	
	cudaMalloc( (void**)&dev_a0, size );
	cudaMalloc( (void**)&dev_b0, size );
	cudaMalloc( (void**)&dev_c0, size );
	cudaMalloc( (void**)&dev_d0, size );
	cudaMalloc( (void**)&dev_A0, size );
	cudaMalloc( (void**)&dev_B0, size );
	cudaMalloc( (void**)&dev_C0, size );
	cudaMalloc( (void**)&dev_D0, size );
	
	cudaMalloc( (void**)&dev_aux1, resultSize );
	cudaMalloc( (void**)&dev_aux2, resultSize );
	cudaMalloc( (void**)&dev_aux_sum1, resultSize );
	cudaMalloc( (void**)&dev_aux_sum2, resultSize );
	cudaMalloc( (void**)&dev_aux_sum3, resultSize );
	cudaMalloc( (void**)&dev_aux_sum4, resultSize );
	
	/*
	Ip_aux = Fp*Hp;
	O resultado é uma matriz 2x2x7.
	
	aA+bC	aB+bD
	cA+dC	cB+dD
	
	*/
	
	// allocate page-locked memory, used to stream
	cudaHostAlloc( (void**)&host_A, size, cudaHostAllocDefault );
	cudaHostAlloc( (void**)&host_B, size, cudaHostAllocDefault );
	cudaHostAlloc( (void**)&host_C, size, cudaHostAllocDefault );
	cudaHostAlloc( (void**)&host_D, size, cudaHostAllocDefault );
	cudaHostAlloc( (void**)&host_a, size, cudaHostAllocDefault );
	cudaHostAlloc( (void**)&host_b, size, cudaHostAllocDefault );
	cudaHostAlloc( (void**)&host_c, size, cudaHostAllocDefault );
	cudaHostAlloc( (void**)&host_d, size, cudaHostAllocDefault );
		
	//host_aux_sum = (float*) malloc(resultSize);
		
	for (int i = 0; i < hlength; i++) {
		host_a[i] = Fp[0][0][i];
		host_b[i] = Fp[0][1][i];
		host_c[i] = Fp[1][0][i];
		host_d[i] = Fp[1][1][i];
		host_A[i] = Hp[0][0][i];
		host_B[i] = Hp[0][1][i];
		host_C[i] = Hp[1][0][i];
		host_D[i] = Hp[1][1][i];
	}
	
/* 00 */
	cudaMemcpyAsync(dev_a0, host_a, size, cudaMemcpyHostToDevice, stream0 );
	cudaMemcpyAsync(dev_b0, host_b, size, cudaMemcpyHostToDevice, stream1 );
	cudaMemcpyAsync(dev_A0, host_A, size, cudaMemcpyHostToDevice, stream0 );
	cudaMemcpyAsync(dev_C0, host_C, size, cudaMemcpyHostToDevice, stream1 );
	
	// clean auxiliary variables
	initializeArray<<<4,4,0,stream0>>>(dev_aux1, resultLength);
	initializeArray<<<4,4,0,stream1>>>(dev_aux2, resultLength);

	// Invoke kernel
    multiplyPolynomialCoefficientsKernel<<<hlength, 1, 0, stream0>>>(dev_a0, hlength, dev_A0, hlength, dev_aux1);
    multiplyPolynomialCoefficientsKernel<<<hlength, 1, 0, stream1>>>(dev_b0, hlength, dev_C0, hlength, dev_aux2);
				
	sumPolynomialCoefficientsKernel<<<resultSize, 1, 0, stream0>>>(dev_aux1, dev_aux2, dev_aux_sum1, resultLength);
	
	
/* 01 */
	// clean auxiliary variables
	initializeArray<<<4,4,0,stream1>>>(dev_aux2, resultLength);
	initializeArray<<<4,4,0,stream0>>>(dev_aux1, resultLength);

	cudaMemcpyAsync(dev_B0, host_B, size, cudaMemcpyHostToDevice, stream0 );
	cudaMemcpyAsync(dev_D0, host_D, size, cudaMemcpyHostToDevice, stream1 );
	
	// Invoke kernel
    multiplyPolynomialCoefficientsKernel<<<hlength, 1, 0, stream1>>>(dev_b0, hlength, dev_D0, hlength, dev_aux2);
    multiplyPolynomialCoefficientsKernel<<<hlength, 1, 0, stream0>>>(dev_a0, hlength, dev_B0, hlength, dev_aux1);
			
	sumPolynomialCoefficientsKernel<<<resultSize, 1, 0, stream1>>>(dev_aux1, dev_aux2, dev_aux_sum2, resultLength);
	

/* 10 */
	// clean auxiliary variables
	initializeArray<<<4,4,0,stream0>>>(dev_aux1, resultLength);
	initializeArray<<<4,4,0,stream1>>>(dev_aux2, resultLength);
	
	cudaMemcpyAsync(dev_c0, host_c, size, cudaMemcpyHostToDevice, stream0 );
	cudaMemcpyAsync(dev_d0, host_d, size, cudaMemcpyHostToDevice, stream1 );
	// Invoke kernel
    multiplyPolynomialCoefficientsKernel<<<hlength, 1, 0, stream0>>>(dev_c0, hlength, dev_A0, hlength, dev_aux1);
    multiplyPolynomialCoefficientsKernel<<<hlength, 1, 0, stream1>>>(dev_d0, hlength, dev_C0, hlength, dev_aux2);
				
	sumPolynomialCoefficientsKernel<<<resultSize, 1, 0, stream0>>>(dev_aux1, dev_aux2, dev_aux_sum3, resultLength);
	
/* 11 */
	// clean auxiliary variables
	initializeArray<<<4,4,0,stream1>>>(dev_aux2, resultLength);
	initializeArray<<<4,4,0,stream0>>>(dev_aux1, resultLength);

	// Invoke kernel
    multiplyPolynomialCoefficientsKernel<<<hlength, 1, 0, stream1>>>(dev_d0, hlength, dev_D0, hlength, dev_aux2);
    multiplyPolynomialCoefficientsKernel<<<hlength, 1, 0, stream0>>>(dev_c0, hlength, dev_B0, hlength, dev_aux1);
			
	sumPolynomialCoefficientsKernel<<<resultSize, 1, 0, stream1>>>(dev_aux1, dev_aux2, dev_aux_sum4, resultLength);
	
	cudaStreamSynchronize( stream0 );
	cudaStreamSynchronize( stream1 );
	
	cudaMemcpy(Ip_aux[0][0], dev_aux_sum1, resultSize, cudaMemcpyDeviceToHost);
	cudaMemcpy(Ip_aux[0][1], dev_aux_sum2, resultSize, cudaMemcpyDeviceToHost);
	cudaMemcpy(Ip_aux[1][0], dev_aux_sum3, resultSize, cudaMemcpyDeviceToHost);
	cudaMemcpy(Ip_aux[1][1], dev_aux_sum4, resultSize, cudaMemcpyDeviceToHost);
	
	cudaEventRecord( stop, 0 );
	cudaEventSynchronize( stop ); 
	cudaEventElapsedTime( &elapsedTime, start, stop ); 
	//printf( "Time taken getIpAux: %3.1f ms\n", elapsedTime );
	
	cudaFreeHost(host_A);
	cudaFreeHost(host_B);
	cudaFreeHost(host_C);
	cudaFreeHost(host_D);
	cudaFreeHost(host_a);
	cudaFreeHost(host_b);
	cudaFreeHost(host_c);
	cudaFreeHost(host_d);
	cudaFree(dev_a0);
	cudaFree(dev_b0);
	cudaFree(dev_c0);
	cudaFree(dev_d0);
	cudaFree(dev_A0);
	cudaFree(dev_B0);
	cudaFree(dev_C0);
	cudaFree(dev_D0);
	cudaStreamDestroy(stream0);
	cudaStreamDestroy(stream1);
	
	return Ip_aux;
}