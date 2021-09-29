#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *plhs = mxGetPr(prhs[0]);
    size_t N = mxGetNumberOfElements(prhs[ 1]);
    double idx_mod = *mxGetPr(prhs[1]);
    double mod_val = *mxGetPr(prhs[2]);
    // You're not supposed to do in-place operations! Don't do this!
    for (ptrdiff_t ii=0; ii<N; ii++) {
        plhs[idx_mod[ii]] += mod_val[ii];
    }
}