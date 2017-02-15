#ifndef __cudaGlobals
#define __cudaGlobals
extern real * devPii[2];
extern real * devTau[2];
extern real * devAReal[2];
extern real * devAImag[2];
extern real * devBReal[2];
extern real * devBImag[2];
extern real * devII[2];
extern int * devNmax[2];
extern cudaStream_t stream[2];
extern cudaStream_t streamRef[2];
extern real * devReferences[2];
extern real * devErr[2];
extern real * devPatterns[2];
extern real * devInvRSquare[2];
extern real * devPSquare[2];
extern real * devMin[2];
extern real * devMax[2];
extern int * devMinIndex[2];
extern real * devMedian[2];

//extern real * devIn[3];
extern int * devOut;
//extern cudaStream_t streamMin[3];
#endif
