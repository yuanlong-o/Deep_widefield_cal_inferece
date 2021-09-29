#include "cuda5s2b.h"
#include "mex.h"
#include <string.h> //for memcpy

#define CUDA5S_LISTDEVICES_MAXDEVICES 256

template <typename T,unsigned S>
        inline unsigned arraysize(const T (&v)[S]) { return S; }

void SetStructFieldToIntRowVector(mxArray *s, const char *fieldname, int *data, int n);
mxArray *makepropstruct(cudaDeviceProp *p);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{   
    if (nrhs != 0 || nlhs != 1) {
        mexErrMsgTxt("no inputs, one output");
    }
    
    cudaDeviceProp pDeviceList[CUDA5S_LISTDEVICES_MAXDEVICES];
    int ndevices = cudaDeviceList_wrapper(pDeviceList, CUDA5S_LISTDEVICES_MAXDEVICES);
    
    if (ndevices < 0) {
        mexPrintf("error code: %d\n", ndevices);
        mexErrMsgTxt("failed to list devices");
    }
    
    mwSize dims[2];
    dims[0] = 1;
    dims[1] = ndevices;
    plhs[0] = mxCreateCellArray(2, dims);
    
    for (int j = 0; j < ndevices; j++) {
        
        mxSetCell(plhs[0], j, makepropstruct(pDeviceList + j));
        
    }
}


mxArray *makepropstruct(cudaDeviceProp *p) {
    mwSize dims[2];
    dims[0] = 1;
    dims[1] = 1;

    const char *field_names[] = { \
            "name", "totalGlobalMem","sharedMemPerBlock","regsPerBlock","warpSize","memPitch","maxThreadsPerBlock","maxThreadsDim", \
            "maxGridSize","clockRate","totalConstMem","major","minor","textureAlignment","texturePitchAlignment","deviceOverlap", \
            "multiProcessorCount","kernelExecTimeoutEnabled", "integrated","canMapHostMemory","computeMode","maxTexture1D", \
            "maxTexture1DMipmap","maxTexture1DLinear","maxTexture2D","maxTexture2DMipmap","maxTexture2DLinear","maxTexture2DGather", \
            "maxTexture3D","maxTexture3DAlt","maxTextureCubemap", "maxTexture1DLayered", "maxTexture2DLayered","maxTextureCubemapLayered", \
            "maxSurface1D","maxSurface2D","maxSurface3D","maxSurface1DLayered","maxSurface2DLayered","maxSurfaceCubemap","maxSurfaceCubemapLayered",\
            "surfaceAlignment","concurrentKernels","ECCEnabled","pciBusID","pciDeviceID","pciDomainID","tccDriver","asyncEngineCount", \
            "unifiedAddressing","memoryClockRate","memoryBusWidth","l2CacheSize","maxThreadsPerMultiProcessor","streamPrioritiesSupported", \
            "globalL1CacheSupported","localL1CacheSupported","sharedMemPerMultiprocessor","regsPerMultiprocessor","managedMemory", "isMultiGpuBoard", "multiGpuBoardGroupID"};

    mxArray *s = mxCreateStructArray(2, dims,
            arraysize(field_names), field_names);    
    
    mxSetField(s, 0, "name", mxCreateString (p->name)); //null terminated.
    mxSetField(s, 0, "maxThreadsPerBlock", mxCreateDoubleScalar (p->maxThreadsPerBlock));
    mxSetField(s, 0, "totalGlobalMem", mxCreateDoubleScalar (p->totalGlobalMem));
    mxSetField(s, 0, "sharedMemPerBlock", mxCreateDoubleScalar (p->sharedMemPerBlock));
    mxSetField(s, 0, "regsPerBlock", mxCreateDoubleScalar (p->regsPerBlock));
    mxSetField(s, 0, "warpSize", mxCreateDoubleScalar (p->warpSize));
    mxSetField(s, 0, "memPitch", mxCreateDoubleScalar (p->memPitch));
    mxSetField(s, 0, "clockRate", mxCreateDoubleScalar (p->clockRate));
    mxSetField(s, 0, "totalConstMem", mxCreateDoubleScalar (p->totalConstMem));
    mxSetField(s, 0, "major", mxCreateDoubleScalar (p->major));
    mxSetField(s, 0, "minor", mxCreateDoubleScalar (p->minor));
    mxSetField(s, 0, "textureAlignment", mxCreateDoubleScalar (p->textureAlignment));
    mxSetField(s, 0, "texturePitchAlignment", mxCreateDoubleScalar (p->texturePitchAlignment));
    mxSetField(s, 0, "deviceOverlap", mxCreateDoubleScalar (p->deviceOverlap));
    mxSetField(s, 0, "multiProcessorCount", mxCreateDoubleScalar (p->multiProcessorCount));
    mxSetField(s, 0, "kernelExecTimeoutEnabled", mxCreateDoubleScalar (p->kernelExecTimeoutEnabled));
    mxSetField(s, 0, "integrated", mxCreateDoubleScalar (p->integrated));
    mxSetField(s, 0, "canMapHostMemory", mxCreateDoubleScalar (p->canMapHostMemory));
    mxSetField(s, 0, "computeMode", mxCreateDoubleScalar (p->computeMode));
    mxSetField(s, 0, "maxTexture1D", mxCreateDoubleScalar (p->maxTexture1D));
    mxSetField(s, 0, "maxTexture1DMipmap", mxCreateDoubleScalar (p->maxTexture1DMipmap));
    mxSetField(s, 0, "maxTexture1DLinear", mxCreateDoubleScalar (p->maxTexture1DLinear));
    mxSetField(s, 0, "maxTextureCubemap", mxCreateDoubleScalar (p->maxTextureCubemap));
    mxSetField(s, 0, "maxSurface1D", mxCreateDoubleScalar (p->maxSurface1D));
    mxSetField(s, 0, "maxSurfaceCubemap", mxCreateDoubleScalar (p->maxSurfaceCubemap));
    mxSetField(s, 0, "surfaceAlignment", mxCreateDoubleScalar (p->surfaceAlignment));
    mxSetField(s, 0, "concurrentKernels", mxCreateDoubleScalar (p->concurrentKernels));
    mxSetField(s, 0, "ECCEnabled", mxCreateDoubleScalar (p->ECCEnabled));
    mxSetField(s, 0, "pciBusID", mxCreateDoubleScalar (p->pciBusID));
    mxSetField(s, 0, "pciDeviceID", mxCreateDoubleScalar (p->pciDeviceID));
    mxSetField(s, 0, "pciDomainID", mxCreateDoubleScalar (p->pciDomainID));
    mxSetField(s, 0, "tccDriver", mxCreateDoubleScalar (p->tccDriver));
    mxSetField(s, 0, "asyncEngineCount", mxCreateDoubleScalar (p->asyncEngineCount));
    mxSetField(s, 0, "unifiedAddressing", mxCreateDoubleScalar (p->unifiedAddressing));
    mxSetField(s, 0, "memoryClockRate", mxCreateDoubleScalar (p->memoryClockRate));
    mxSetField(s, 0, "memoryBusWidth", mxCreateDoubleScalar (p->memoryBusWidth));
    mxSetField(s, 0, "l2CacheSize", mxCreateDoubleScalar (p->l2CacheSize));
    mxSetField(s, 0, "maxThreadsPerMultiProcessor", mxCreateDoubleScalar (p->maxThreadsPerMultiProcessor));
    mxSetField(s, 0, "streamPrioritiesSupported", mxCreateDoubleScalar (p->streamPrioritiesSupported));
    mxSetField(s, 0, "globalL1CacheSupported", mxCreateDoubleScalar (p->globalL1CacheSupported));
    mxSetField(s, 0, "localL1CacheSupported", mxCreateDoubleScalar (p->localL1CacheSupported));
    mxSetField(s, 0, "sharedMemPerMultiprocessor", mxCreateDoubleScalar (p->sharedMemPerMultiprocessor));
    mxSetField(s, 0, "regsPerMultiprocessor", mxCreateDoubleScalar (p->regsPerMultiprocessor));
    mxSetField(s, 0, "managedMemory", mxCreateDoubleScalar (p->managedMemory));
    mxSetField(s, 0, "isMultiGpuBoard", mxCreateDoubleScalar (p->isMultiGpuBoard));
    mxSetField(s, 0, "multiGpuBoardGroupID", mxCreateDoubleScalar (p->multiGpuBoardGroupID));
    
    SetStructFieldToIntRowVector(s, "maxTexture2D", p->maxTexture2D, arraysize(p->maxTexture2D));
    SetStructFieldToIntRowVector(s, "maxTexture3D", p->maxTexture2D, arraysize(p->maxTexture3D));
    SetStructFieldToIntRowVector(s, "maxSurfaceCubemapLayered", p->maxSurfaceCubemapLayered, arraysize(p->maxSurfaceCubemapLayered));
    SetStructFieldToIntRowVector(s, "maxTexture2DMipmap", p->maxTexture2DMipmap, arraysize(p->maxTexture2DMipmap));
    SetStructFieldToIntRowVector(s, "maxTexture2DLinear", p->maxTexture2DLinear, arraysize(p->maxTexture2DLinear));
    SetStructFieldToIntRowVector(s, "maxTexture2DGather", p->maxTexture2DGather, arraysize(p->maxTexture2DGather));
    SetStructFieldToIntRowVector(s, "maxTexture3DAlt", p->maxTexture3DAlt, arraysize(p->maxTexture3DAlt));
    SetStructFieldToIntRowVector(s, "maxTexture1DLayered", p->maxTexture1DLayered, arraysize(p->maxTexture1DLayered));
    SetStructFieldToIntRowVector(s, "maxTexture2DLayered", p->maxTexture2DLayered, arraysize(p->maxTexture2DLayered));
    SetStructFieldToIntRowVector(s, "maxTextureCubemapLayered", p->maxTextureCubemapLayered, arraysize(p->maxTextureCubemapLayered));
    SetStructFieldToIntRowVector(s, "maxSurface2D", p->maxSurface2D, arraysize(p->maxSurface2D));
    SetStructFieldToIntRowVector(s, "maxSurface3D", p->maxSurface3D, arraysize(p->maxSurface3D));
    SetStructFieldToIntRowVector(s, "maxSurface1DLayered", p->maxSurface1DLayered, arraysize(p->maxSurface1DLayered));
    SetStructFieldToIntRowVector(s, "maxSurface2DLayered", p->maxSurface2DLayered, arraysize(p->maxSurface2DLayered));
    SetStructFieldToIntRowVector(s, "maxThreadsDim", p->maxThreadsDim, arraysize(p->maxThreadsDim));
    SetStructFieldToIntRowVector(s, "maxGridSize", p->maxGridSize, arraysize(p->maxGridSize));
    
    return s;
}

void SetStructFieldToIntRowVector(mxArray *s, const char *fieldname, int *data, int n) {
    mxArray *f = mxCreateNumericMatrix(1, n, mxINT32_CLASS, mxREAL);
    memcpy(mxGetData(f), data, sizeof(int) * n);
    mxSetField(s, 0, fieldname, f);
}
