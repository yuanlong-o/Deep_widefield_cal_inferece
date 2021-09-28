//usage: cuda5s_mex_clear(g)
#include "cuda5s2b.h"
#include "mex.h"
#include "matrix.h"
#include <string.h>

#define VAR_gdata prhs[0]

#define NINPUTS 1
#define NOUTPUTS 0

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
    if (nrhs != NINPUTS || nlhs != NOUTPUTS) {
        mexErrMsgTxt("must have NINPUTS input(s) and NOUTPUTS output(s)");
    }
    if (mxUINT8_CLASS != mxGetClassID(VAR_gdata) || mxGetNumberOfElements(VAR_gdata) != sizeof(gpu5s_problem)) {
        mexErrMsgTxt("g must be a uint8 array with sizeof(gpu5s_problem) elements");
    }
    gpu5s_problem * gin = (gpu5s_problem *) mxGetData(VAR_gdata);
    free_gpupfdata(gin); //frees GPU arrays and sets gpu array pointers to NULL. does not free gin or any hosts-side arrays
}
