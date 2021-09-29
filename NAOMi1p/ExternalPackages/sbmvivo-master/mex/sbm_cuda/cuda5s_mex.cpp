//usage: [marglik, log_sum_raw_w, neff, nmean] = cuda5s_mex(g, f, u, q_spike, seedval)
#include "cuda5s2b.h"
#include "mex.h"
#include "matrix.h"
#include <cmath>

#define VAR_marglik plhs[0]
#define VAR_log_sum_raw_w plhs[1]
#define VAR_neff     plhs[2]
#define VAR_nmean    plhs[3]
#define VAR_gpmean   plhs[4]
#define VAR_gpsqmean plhs[5]

#define VAR_gdata prhs[0]
#define VAR_f prhs[1]
#define VAR_u prhs[2]
#define VAR_q_spike prhs[3]
#define VAR_rngseed prhs[4]

#define NINPUTS_MIN 4
#define NINPUTS_MAX 5
#define NOUTPUTS 6

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
    if (nrhs < NINPUTS_MIN || nrhs > NINPUTS_MAX || nlhs > NOUTPUTS) {
        mexErrMsgTxt("must have NINPUTS input(s) and can have up to NOUTPUTS output(s)");
    }
    if (mxUINT8_CLASS != mxGetClassID(VAR_gdata) || mxGetNumberOfElements(VAR_gdata) != sizeof(gpu5s_problem)) {
        mexErrMsgTxt("g must be a uint8 array with sizeof(gpu5s_problem) elements");
    }
    gpu5s_problem * g = (gpu5s_problem *) mxGetData(VAR_gdata);
    if (mxSINGLE_CLASS != mxGetClassID(VAR_f) || mxSINGLE_CLASS != mxGetClassID(VAR_u) || mxSINGLE_CLASS != mxGetClassID(VAR_q_spike)) {
        mexErrMsgTxt("f, u and q_spike must be single");
    }
    if (mxGetNumberOfElements(VAR_f) != g->T || mxGetNumberOfElements(VAR_u) != g->T) {
        mexErrMsgTxt("f and u must have T elemnents");
    }
    g->h.fobs = (float *) mxGetData(VAR_f);
    g->h.u = (float *) mxGetData(VAR_u);
    int nq = (g->T + g->options.ntimepoints_pre) * g->options.nsteps;
    if (mxGetNumberOfElements(VAR_q_spike) == 1)
        fillarray(g->d.q_spike, (float) mxGetScalar(VAR_q_spike), nq);
    else {
        if (mxGetNumberOfElements(VAR_q_spike) == nq)
            cudaMemcpy_h2d_wrapper((void *) g->d.q_spike, mxGetData(VAR_q_spike), nq * sizeof(float));
        else
            mexErrMsgTxt("q_spike must have either 1 element or (T + ntimepoints_pre) * nsteps elements");
    }
    if (nrhs > 4) {
        if (mxGetNumberOfElements(VAR_rngseed) != 1 || mxUINT64_CLASS != mxGetClassID(VAR_rngseed))
            mexErrMsgTxt("seedval must be a uint64 scalar");
        unsigned long long seedval = *((unsigned long long *) (mxGetData(VAR_rngseed)));
        int nrng = g->options.nblocks_pf * BLOCKSIZE_PF;
        cudaseedrng(g->options.nblocks_pf, g->d.rngstates, nrng, seedval);
    }
    g->options.computenmean = nlhs > 3;
    pushparameterstodevice(g);
    gpu5s_initialstates(g);
    //mexPrintf("\n*** starting main loop ***\n");
    int mainloopresult = runpfmainloop(g);
    g->h.fobs = NULL; g->h.u = NULL;
    if (mainloopresult < 0)
        mexErrMsgTxt("main loop failed");
    if (!nlhs) return;
    if (nlhs == 1) {
        VAR_marglik = mxCreateDoubleScalar((double) gpu5s_marglik(g));
        return;
    }
    
    VAR_log_sum_raw_w = mxCreateNumericMatrix(1, g->T, mxSINGLE_CLASS, mxREAL);
    float * plsrw = (float *) mxGetData(VAR_log_sum_raw_w);
    cudaMemcpy_d2h_wrapper((void *) plsrw, (void *) g->d.log_sum_raw_w, g->T * sizeof(float));
    float logmarglik = 0.0;
    for (int t = 0; t < g->T; t++)
        logmarglik += plsrw[t];
    VAR_marglik = mxCreateDoubleScalar((double) logmarglik);
    if (nlhs == 2) return;
    
    VAR_neff = mxCreateNumericMatrix(1, g->T, mxSINGLE_CLASS, mxREAL);
    cudaMemcpy_d2h_wrapper(mxGetData(VAR_neff), (void *) g->d.neff, g->T * sizeof(float));
    if (nlhs == 3) return;
    
    //initialize moment outputs
    VAR_nmean    = mxCreateNumericMatrix(1,  nq, mxSINGLE_CLASS, mxREAL);
    VAR_gpmean   = mxCreateNumericMatrix(1,  g->T, mxSINGLE_CLASS, mxREAL);
    VAR_gpsqmean = mxCreateNumericMatrix(1,  g->T, mxSINGLE_CLASS, mxREAL);
    
    //copy moments from GPU arrays to Matlab variables on the CPU
    cudaMemcpy_d2h_wrapper(mxGetData(VAR_nmean),    (void *) g->d.nmean,    nq * sizeof(float));
    cudaMemcpy_d2h_wrapper(mxGetData(VAR_gpmean),   (void *) g->d.gpmean,   g->T * sizeof(float));
    cudaMemcpy_d2h_wrapper(mxGetData(VAR_gpsqmean), (void *) g->d.gpsqmean, g->T * sizeof(float));
    
}
