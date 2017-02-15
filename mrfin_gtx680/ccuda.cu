#include"globals.h"
#include"cudaGlobals.h"
#include"cudaError.h"
real * devPii[2];
real * devTau[2];
real * devAReal[2];
real * devAImag[2];
real * devBReal[2];
real * devBImag[2];
real * devII[2];
int * devNmax[2];
cudaStream_t stream[2];
real * devReferences[2];
real * devErr[2];
real * devPatterns[2];
real * devInvRSquare[2];
real * devPSquare[2];
cudaStream_t streamRef[2];
real * devMin[2];
real * devMax[2];
int * devMinIndex[2];
real * devMedian[2];
int * devOut;


void freeCudaPointer(void ** pointer) {
	CudaSafeCall(cudaFreeHost(*pointer));
}
void allocCudaPointer(void ** pointer, size_t size) {
	CudaSafeCall(cudaMallocHost((void**)pointer, size));
}
	
void mallocCudaReferences(int i, int const mPatterns, int const nPatterns, int const mReferences, int const nReferences ) {
			CudaSafeCall(cudaStreamCreate(&streamRef[i]));
			CudaSafeCall(cudaMalloc((void**)&devPatterns[i], mPatterns*nPatterns*sizeof(real)));
			CudaSafeCall(cudaMalloc((void**)&devReferences[i],mReferences*nReferences*sizeof(real)));
			CudaSafeCall(cudaMalloc((void**)&devInvRSquare[i], mReferences*sizeof(real)));
			CudaSafeCall(cudaMalloc((void**)&devPSquare[i], mPatterns*sizeof(real)));
			CudaSafeCall(cudaMalloc((void**)&devErr[i], mPatterns*mReferences*sizeof(real)));
			CudaSafeCall(cudaMalloc((void**)&devMin[i], mPatterns*sizeof(real)));
			CudaSafeCall(cudaMalloc((void**)&devMax[i], mPatterns*sizeof(real))); //TODO: czy rozmiar dobry? (04.04.13 by szmigacz)
			CudaSafeCall(cudaMalloc((void**)&devMinIndex[i], mPatterns*sizeof(int)));
			CudaSafeCall(cudaMalloc((void**)&devMedian[i], mPatterns*sizeof(real)));

}

void freeCudaMemory() {
	#ifdef CUDA
		for(int i=0;i<2;i++) {
			CudaSafeCall(cudaStreamSynchronize(stream[i]));
		}
	
		for(int i=0;i<2;i++) {
			CudaSafeCall(cudaFree(devPii[i]));
			CudaSafeCall(cudaFree(devTau[i]));
			CudaSafeCall(cudaFree(devAReal[i]));
			CudaSafeCall(cudaFree(devAImag[i]));
			CudaSafeCall(cudaFree(devBReal[i]));
			CudaSafeCall(cudaFree(devBImag[i]));
			CudaSafeCall(cudaFree(devII[i]));
			CudaSafeCall(cudaFree(devNmax[i]));
			CudaSafeCall(cudaStreamDestroy(stream[i]));
		}
	#endif //CUDA
}

void freeCudaRefMemory() {
	#ifdef CUDA
		for(int i=0;i<2;i++) {
			CudaSafeCall(cudaStreamSynchronize(streamRef[i]));
		}
		for(int i=0;i<2;i++) {
	
			CudaSafeCall(cudaFree(devInvRSquare[i]));
			CudaSafeCall(cudaFree(devPSquare[i]));
			CudaSafeCall(cudaFree(devErr[i]));
			CudaSafeCall(cudaFree(devPatterns[i]));
			CudaSafeCall(cudaFree(devReferences[i]));
			CudaSafeCall(cudaFree(devMin[i]));
			CudaSafeCall(cudaFree(devMax[i]));
			CudaSafeCall(cudaFree(devMinIndex[i]));
			CudaSafeCall(cudaFree(devMedian[i]));
			CudaSafeCall(cudaStreamDestroy(streamRef[i]));
		}
	#endif //CUDA
}

void cudaFinalize() {
	#ifdef CUDA
		CudaSafeCall(cudaDeviceSynchronize());
	#endif //CUDA
}
void cuda1stPolarizationSync() {
	#ifdef CUDA
		CudaSafeCall(cudaStreamSynchronize(streamRef[0]));
		CudaSafeCall(cudaStreamSynchronize(streamRef[1]));
	#endif //CUDA
}

void freeCudaMemoryMin() {
	#ifdef CUDA
			CudaSafeCall(cudaFree(devOut));
	#endif //CUDA
}
