#include<cstdio>
#include<omp.h>
#include"globals.h"
#include"cudaGlobals.h"
#include"cudaError.h"
#define sq(x) ((x)*(x))
//tutaj tez zrobic te optymalizacje zeby 2x nie czytac najlepiej wzgledem wymiaru Xowego
#define THREADS 128


texture<float> texPii;
texture<float> texTau;
texture<float> texaReal;
texture<float> texbReal;
texture<float> texaImag;
texture<float> texbImag;

#if __CUDA_ARCH__ >= 300
__global__ void kernelGenerate(int const rSize, int const nPiiTau, real const * const pii, real const * const tau, 
		real const * const aReal, real const * const aImag, real const * const bReal, real const * const bImag,
		int const * const NmaxTable, real * const II) {
	volatile __shared__ real realis[4];
	volatile __shared__ real imaginalis[4];
	const int Nmax = NmaxTable[blockIdx.y];
	int index  = blockIdx.x * nPiiTau + threadIdx.x;
	int indexY = blockIdx.y * nPiiTau + threadIdx.x;
	real r=0.0f;
	real i=0.0f;
	if(threadIdx.x < Nmax) {
		//const real pi = pii[index];
		//const real ta = tau[index];
		const real pi = tex1Dfetch(texPii, index);
		const real ta = tex1Dfetch(texTau, index);
		//r = aReal[indexY] * pi + bReal[indexY] * ta;
		//i = aImag[indexY] * pi + bImag[indexY] * ta;
		r = tex1Dfetch(texaReal,indexY) * pi + tex1Dfetch(texbReal,indexY) * ta;
		i = tex1Dfetch(texaImag,indexY) * pi + tex1Dfetch(texbImag,indexY) * ta;
		for(int id = THREADS ; id < Nmax ; id+=THREADS) {
			if(threadIdx.x + id < Nmax) {
				index += THREADS; 
				indexY += THREADS;
				//const real pi = pii[index];
				//const real ta = tau[index];
				const real pi = tex1Dfetch(texPii, index);
				const real ta = tex1Dfetch(texTau, index);
				//r += aReal[indexY] * pi + bReal[indexY] * ta;
				//i += aImag[indexY] * pi + bImag[indexY] * ta;
				r += tex1Dfetch(texaReal,indexY) * pi + tex1Dfetch(texbReal,indexY) * ta;
				i += tex1Dfetch(texaImag,indexY) * pi + tex1Dfetch(texbImag,indexY) * ta;
			}
		}
	}

	//butterfly reduction across warp
	for (int j=16; j>=1; j/=2) {
		r += __shfl_xor(r , j, 32);
		i += __shfl_xor(i , j, 32);
	}
	//further reduction across block
	if(threadIdx.x % 32 == 0) {
		realis[threadIdx.x>>5] = r;
		imaginalis[threadIdx.x>>5] = i;
	}
	__syncthreads();
	if(threadIdx.x <2) {
		realis[threadIdx.x] += realis[threadIdx.x+2];
		imaginalis[threadIdx.x] += imaginalis[threadIdx.x+2];
		realis[threadIdx.x] += realis[threadIdx.x+1];
		imaginalis[threadIdx.x] += imaginalis[threadIdx.x+1];
	}

	if(threadIdx.x==0)
		II[blockIdx.x + blockIdx.y*gridDim.x] = sq(realis[0]) + sq(imaginalis[0]);
}
#endif
#if __CUDA_ARCH__ < 300
__global__ void kernelGenerate(int const rSize, int const nPiiTau, real const * const pii, real const * const tau, real const * const aReal, real const * const aImag, real const * const bReal, real const * const bImag, int const * const NmaxTable, real * const II) {
	volatile __shared__ real realis[THREADS];
	volatile __shared__ real imaginalis[THREADS];
	const int Nmax = NmaxTable[blockIdx.y];
	int index = blockIdx.x * nPiiTau + threadIdx.x;
	int indexY = blockIdx.y * nPiiTau + threadIdx.x; 
	real r=0.0f;
	real i=0.0f;
	if(threadIdx.x < Nmax) {
		const real pi = pii[index];
		const real ta = tau[index];
		//const real pi = tex1Dfetch(texPii, index);
		//const real ta = tex1Dfetch(texTau, index);
		r = aReal[indexY] * pi + bReal[indexY] * ta;
		i = aImag[indexY] * pi + bImag[indexY] * ta;
		//r = tex1Dfetch(texaReal,indexY) * pi + tex1Dfetch(texbReal,indexY) * ta;
		//i = tex1Dfetch(texaImag,indexY) * pi + tex1Dfetch(texbImag,indexY) * ta;
		for(int id = THREADS ; id < Nmax ; id+=THREADS) {
			if(threadIdx.x + id < Nmax) {
				index += THREADS; 
				indexY += THREADS;
				const real pi = pii[index];
				const real ta = tau[index];
				//const real pi = tex1Dfetch(texPii, index);
				//const real ta = tex1Dfetch(texTau, index);
				r += aReal[indexY] * pi + bReal[indexY] * ta;
				i += aImag[indexY] * pi + bImag[indexY] * ta;
				//r += tex1Dfetch(texaReal,indexY) * pi + tex1Dfetch(texbReal,indexY) * ta;
				//i += tex1Dfetch(texaImag,indexY) * pi + tex1Dfetch(texbImag,indexY) * ta;
			}
		}
	}
	if ( threadIdx.x < Nmax ) {
		realis[threadIdx.x] = r;
		imaginalis[threadIdx.x] = i;
	}
	else {
		realis[threadIdx.x]=0.0f;
		imaginalis[threadIdx.x]=0.0f;
	}
	__syncthreads();
	if(threadIdx.x < 64 ) {
		realis[threadIdx.x]+=realis[threadIdx.x+64];
		imaginalis[threadIdx.x]+=imaginalis[threadIdx.x+64];
	}
	__syncthreads();
	if(threadIdx.x < 32) {
		realis[threadIdx.x]+=realis[threadIdx.x+32];
		realis[threadIdx.x]+=realis[threadIdx.x+16];
		realis[threadIdx.x]+=realis[threadIdx.x+8];
		realis[threadIdx.x]+=realis[threadIdx.x+4];
		realis[threadIdx.x]+=realis[threadIdx.x+2];
		realis[threadIdx.x]+=realis[threadIdx.x+1];
		imaginalis[threadIdx.x]+=imaginalis[threadIdx.x+32];
		imaginalis[threadIdx.x]+=imaginalis[threadIdx.x+16];
		imaginalis[threadIdx.x]+=imaginalis[threadIdx.x+8];
		imaginalis[threadIdx.x]+=imaginalis[threadIdx.x+4];
		imaginalis[threadIdx.x]+=imaginalis[threadIdx.x+2];
		imaginalis[threadIdx.x]+=imaginalis[threadIdx.x+1];
		if(threadIdx.x==0)
			II[blockIdx.x + blockIdx.y*gridDim.x] = sq(realis[0]) + sq(imaginalis[0]);
	}
}


#endif

void cudaGenerate(int rSize, int pattern_length, int * Nmax, real * pii, int nPiiTau,  real * tau, real * aReal, real * aImag, real * bReal, real * bImag, real * II, int polarization ) {
	CudaSafeCall(cudaStreamCreate(&stream[polarization]));
	CudaSafeCall(cudaMalloc((void**)&devPii[polarization]  , nPiiTau*pattern_length*sizeof(real)));
	CudaSafeCall(cudaMalloc((void**)&devTau[polarization]  , nPiiTau*pattern_length*sizeof(real)));
	CudaSafeCall(cudaMalloc((void**)&devAReal[polarization], rSize*nPiiTau*sizeof(real)));
	CudaSafeCall(cudaMalloc((void**)&devBReal[polarization], rSize*nPiiTau*sizeof(real)));
	CudaSafeCall(cudaMalloc((void**)&devAImag[polarization], rSize*nPiiTau*sizeof(real)));
	CudaSafeCall(cudaMalloc((void**)&devBImag[polarization], rSize*nPiiTau*sizeof(real)));
	CudaSafeCall(cudaMalloc((void**)&devII[polarization]   , rSize*pattern_length*sizeof(real)));
	CudaSafeCall(cudaMalloc((void**)&devNmax[polarization] , rSize*sizeof(int)));

	CudaSafeCall(cudaMemcpyAsync(devPii[polarization]  , pii  , nPiiTau*pattern_length*sizeof(real) , cudaMemcpyHostToDevice , stream[polarization]));
	CudaSafeCall(cudaMemcpyAsync(devTau[polarization]  , tau  , nPiiTau*pattern_length*sizeof(real) , cudaMemcpyHostToDevice , stream[polarization]));
	CudaSafeCall(cudaMemcpyAsync(devNmax[polarization] , Nmax , rSize*sizeof(int)                   , cudaMemcpyHostToDevice , stream[polarization]));
	CudaSafeCall( cudaBindTexture( NULL, texPii, devPii[polarization], nPiiTau*pattern_length*sizeof(real)));
	CudaSafeCall( cudaBindTexture( NULL, texTau, devTau[polarization], nPiiTau*pattern_length*sizeof(real)));

	CudaSafeCall(cudaMemcpyAsync(devAReal[polarization], aReal, rSize*nPiiTau*sizeof(real), cudaMemcpyHostToDevice, stream[polarization]));
	CudaSafeCall(cudaMemcpyAsync(devBReal[polarization], bReal, rSize*nPiiTau*sizeof(real), cudaMemcpyHostToDevice, stream[polarization]));
	CudaSafeCall(cudaMemcpyAsync(devAImag[polarization], aImag, rSize*nPiiTau*sizeof(real), cudaMemcpyHostToDevice, stream[polarization]));
	CudaSafeCall(cudaMemcpyAsync(devBImag[polarization], bImag, rSize*nPiiTau*sizeof(real), cudaMemcpyHostToDevice, stream[polarization]));
	CudaSafeCall( cudaBindTexture( NULL, texaReal, devAReal[polarization], rSize*nPiiTau*sizeof(real)));
	CudaSafeCall( cudaBindTexture( NULL, texbReal, devBReal[polarization], rSize*nPiiTau*sizeof(real)));
	CudaSafeCall( cudaBindTexture( NULL, texaImag, devAImag[polarization], rSize*nPiiTau*sizeof(real)));
	CudaSafeCall( cudaBindTexture( NULL, texbImag, devBImag[polarization], rSize*nPiiTau*sizeof(real)));

#ifdef GF580
	cudaFuncSetCacheConfig(kernelGenerate, cudaFuncCachePreferL1);
#endif //GF580
#ifdef GF680
	cudaFuncSetCacheConfig(kernelGenerate, cudaFuncCachePreferShared);
#endif //GF680
	kernelGenerate<<<dim3(pattern_length,rSize,1), THREADS, 0, stream[polarization]>>>(rSize, nPiiTau, devPii[polarization], devTau[polarization],
			devAReal[polarization], devAImag[polarization], devBReal[polarization], devBImag[polarization],
			devNmax[polarization], devII[polarization]);

	CudaSafeCall(cudaMemcpyAsync(II, devII[polarization], rSize*pattern_length*sizeof(real), cudaMemcpyDeviceToHost, stream[polarization]));
	CudaSafeCall(cudaUnbindTexture( texPii));
	CudaSafeCall(cudaUnbindTexture( texTau));
	CudaSafeCall(cudaUnbindTexture( texaReal));
	CudaSafeCall(cudaUnbindTexture( texbReal));
	CudaSafeCall(cudaUnbindTexture( texbImag));
	CudaSafeCall(cudaUnbindTexture( texaImag));
}
#undef sq



