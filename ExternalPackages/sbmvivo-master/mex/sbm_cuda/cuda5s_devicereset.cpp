#include "cuda5s2b.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
    cudaDeviceReset_wrapper();
}
