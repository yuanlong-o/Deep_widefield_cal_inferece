//usage: g = cuda5s_mex_configure(g, np, P, dt, nA2D, T, nsteps, ntimepoints_pre, K, resamplethreshold, n_newton_iterations, FBGp1)
#include "cuda5s2b.h"
#include "mex.h"
#include "matrix.h"
#include <string.h>
#include <math.h>

#define VAR_goutput plhs[0]
#define VAR_ginput prhs[0]
#define VAR_nparticles prhs[1]
#define VAR_P prhs[2]
#define VAR_dt prhs[3]
#define VAR_nA2D prhs[4]
#define VAR_T prhs[5]
#define VAR_nsteps prhs[6]
#define VAR_ntimepoints_pre prhs[7]
#define VAR_K prhs[8]
#define VAR_resamplethreshold prhs[9]
#define VAR_n_newton_iterations prhs[10]
#define VAR_FBGp1 prhs[11]
#define NINPUTS 12
#define NOUTPUTS 1

float pval(const mxArray *s, const char *fieldname);
float pvalel(const mxArray *s, const char *fieldname, unsigned int ii);

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
    if (nrhs != NINPUTS || nlhs != NOUTPUTS)
        mexErrMsgTxt("must have NINPUTS input(s) and NOUTPUTS output(s)");
    if (!mxIsNumeric(VAR_dt) || !mxIsNumeric(VAR_nA2D) || !mxIsNumeric(VAR_resamplethreshold) || !mxIsNumeric(VAR_nparticles))
        mexErrMsgTxt("dt, nA2D, nsteps and resamplethreshold must be numeric");
    if (mxSTRUCT_CLASS != mxGetClassID(VAR_P))
        mexErrMsgTxt("P must be a structure");
    if (mxGetNumberOfElements(VAR_dt) != 1 || mxGetNumberOfElements(VAR_P) != 1 || mxGetNumberOfElements(VAR_nA2D) != 1)
        mexErrMsgTxt("dt, P, V, nA2D and resamplethreshold must have one element");
    if (mxUINT8_CLASS != mxGetClassID(VAR_ginput) || mxGetNumberOfElements(VAR_ginput) != sizeof(gpu5s_problem))
        mexErrMsgTxt("g must be a uint8 array with sizeof(gpu5s_problem) elements");
    float nA2D = (float) mxGetScalar(VAR_nA2D);
    double np_d = mxGetScalar(VAR_nparticles);
    int np = (int) np_d;
    int T = (int) mxGetScalar(VAR_T);
    int K = (int) mxGetScalar(VAR_K);
    int nsteps = (int) mxGetScalar(VAR_nsteps);
    int ntimepoints_pre = (int) mxGetScalar(VAR_ntimepoints_pre);
    int n_newton_iterations = (int) mxGetScalar(VAR_n_newton_iterations);
    double FBGp1 = mxGetScalar(VAR_FBGp1);
    
    gpu5s_problem * ginput = (gpu5s_problem *) mxGetData(VAR_ginput);
    if (    K                               > ginput->d.maxK             ||
            (K + 1) * nsteps                > ginput->d.maxKp1substeps   ||
            (T + ntimepoints_pre) * nsteps  > ginput->d.maxtotalsubsteps ||
            nsteps                          > ginput->d.maxnsteps        ||
            nsteps * ntimepoints_pre        > ginput->d.maxpresubsteps   ||
            T                               > ginput->d.maxT )
        mexErrMsgTxt("Insufficient space has been allocated for these settings");
    if (((double) (MAX(ntimepoints_pre, K + 1) * nsteps)) * np_d  > 2147483646.5)
        mexErrMsgTxt("Maximum allowed valued of max(ntimepoints_pre, K+1) * nsteps * Nparticles is 2^31 - 1 = 2147483647");
    VAR_goutput = mxCreateNumericMatrix(1, sizeof(gpu5s_problem), mxUINT8_CLASS, mxREAL);
    gpu5s_problem * g = (gpu5s_problem *) mxGetData(VAR_goutput);
    memcpy((void *) g, (void *) ginput, sizeof(gpu5s_problem));
    
    g->dt = (float) mxGetScalar(VAR_dt);
    g->options.resamplethreshold = (float) mxGetScalar(VAR_resamplethreshold);
    g->T = T;
    g->options.K = K;
    g->options.nsteps = nsteps;
    g->options.ntimepoints_pre = ntimepoints_pre;
    g->options.n_newton_iterations = n_newton_iterations;
    if (setgridsize(g, np) < 0) //assigns g->options.nparticles as well as nblocks, nblocks_pf and nblocks_rs
        mexErrMsgTxt("Insufficient space has been allocated for this particle count");
        
    float A = pval(VAR_P, "A");
    g->params.vF = pval(VAR_P, "zeta") / nA2D;
    g->params.sigma_r = pval(VAR_P, "sigma_r");
    g->params.lambda = pval(VAR_P, "fr") * g->dt;
    g->params.S = pval(VAR_P, "S") / A;
    g->params.gain = pval(VAR_P, "gain") / nA2D;
    if (USEKDEX) {
    	g->params.kd_ex = pval(VAR_P, "kd_ex") / A;
    	g->params.maxex = g->params.kd_ex / pval(VAR_P, "tau_ex");
    } else {
    	g->params.maxex = 1.0 / pval(VAR_P, "tau_ex");
    }
    g->params.c0 = pval(VAR_P, "c0") / A;    
    g->params.fdc = pval(VAR_P, "fdc");
    
    g->params.db1 = pvalel(VAR_P, "dbrightness", 0);
    g->params.db2 = pvalel(VAR_P, "dbrightness", 1);
    g->params.db3 = pvalel(VAR_P, "dbrightness", 2);
    g->params.db4 = pvalel(VAR_P, "dbrightness", 3);
    
    g->params.kon0 = pvalel(VAR_P, "kon", 0) * A;
    g->params.kon1 = pvalel(VAR_P, "kon", 1) * A;
    g->params.kon2 = pvalel(VAR_P, "kon", 2) * A;
    g->params.kon3 = pvalel(VAR_P, "kon", 3) * A;
    
    g->params.koff0 = pvalel(VAR_P, "koff", 0);
    g->params.koff1 = pvalel(VAR_P, "koff", 1);
    g->params.koff2 = pvalel(VAR_P, "koff", 2);
    g->params.koff3 = pvalel(VAR_P, "koff", 3);
    
    g->params.kon_B0 = 1.0 / ((pvalel(VAR_P, "kd_B", 0) / A) * pvalel(VAR_P, "tau_B", 0));
    g->params.kon_B1 = 1.0 / ((pvalel(VAR_P, "kd_B", 1) / A) * pvalel(VAR_P, "tau_B", 1));
    g->params.koff_B0 = 1.0 / pvalel(VAR_P, "tau_B", 0);
    g->params.koff_B1 = 1.0 / pvalel(VAR_P, "tau_B", 1);
    
    g->params.FBGp1 = FBGp1;
    
    g->params.Btot0 = pvalel(VAR_P, "Btot", 0) / A;
    g->params.Btot1 = pvalel(VAR_P, "Btot", 1) / A;

    //so far we've just allocated and constructed the gpu5s_problem structure and returned a pointer to it.
    //the parameters will be pushed to the GPU itself when we call cuda5s_mex
}

float pvalel(const mxArray *s, const char *fieldname, unsigned int ii) {
    mxArray* pfield = mxGetField(s, 0, fieldname);
    if (pfield == NULL) {
        mexPrintf(" %s :", fieldname);
        mexErrMsgTxt("Invalid field");
    }
    void * pData = mxGetData(pfield);
    if (mxSINGLE_CLASS == mxGetClassID(pfield)) {
        return (float) *((float *) pData + ii);
    }
    if (mxDOUBLE_CLASS == mxGetClassID(pfield)) {
        return (float) *((double *) pData + ii);
    }
    mexErrMsgTxt("Field must be an array of type single or double");
    return -1.0;
}

float pval(const mxArray *s, const char *fieldname) {
    mxArray* pfield = mxGetField(s, 0, fieldname);
    if (pfield == NULL) {
        mexPrintf(" %s :", fieldname);
        mexErrMsgTxt("Invalid field");
    }
    if (!mxIsNumeric(pfield)) {
        mexErrMsgTxt("Field value must be a numeric array");
    }
    return (float) mxGetScalar(pfield);
}
