#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    //if (!mxIsDouble(prhs[0]) && !mxIsSingle(prhs[0])){
    if (!mxIsSingle(prhs[0])){    
    //    mexErrMsgTxt("Array must be double or single.");
        mexErrMsgTxt("Array must be single.");
    }
    
    if (!mxIsInt32(prhs[1])){
        mexErrMsgTxt("Array must be int32.");
    }
    if (!mxIsSingle(prhs[2])){
        mexErrMsgTxt("Array must be single.");
    }
    //if (mxIsDouble(prhs[0])){
    //    double *A = mxGetPr(prhs[0]);
    //} else {
    //    float *A = (float*)mxGetPr(prhs[0]);
    //}
    //if (mxIsDouble(prhs[1])){
    //    double *idx_mod = mxGetPr(prhs[1]);
    //} else {
    //    float *idx_mod = (float*)mxGetPr(prhs[1]);
    //}
    //if (mxIsDouble(prhs[2])){
    //    double *mod_val = mxGetPr(prhs[2]);
    //} else {
    //    float *mod_val = (float*)mxGetPr(prhs[0]);
    //}
    float *A = (float*)mxGetPr(prhs[0]);
    int32_T *idx_mod = (int32_T*)mxGetPr(prhs[1]);
    float *mod_val = (float*)mxGetPr(prhs[2]);
    float sc = (float)mxGetScalar(prhs[3]);
    long int ii = 0;
    long int kk = 0;
    size_t N = mxGetNumberOfElements(prhs[1]);
    for (ii=0; ii<N; ii++) {
        A[idx_mod[ii]-1] = mod_val[ii]*sc;
    }
}