//usage: g = cuda5s_mex_init(maxT, np, maxnsteps, maxK, maxKp1substeps, maxtotalsubsteps, maxpresubsteps, seedval)
#include "cuda5s2b.h"
#include "mex.h"
#include "matrix.h"
#include <string>

#define VAR_gdata plhs[0]

#define VAR_maxT prhs[0]
#define VAR_np prhs[1]
#define VAR_maxnsteps prhs[2]
#define VAR_maxK prhs[3]
#define VAR_maxKp1substeps prhs[4]
#define VAR_maxtotalsubsteps prhs[5]
#define VAR_maxpresubsteps prhs[6]
#define VAR_seedval prhs[7]
#define NINPUTS 8
#define NOUTPUTS 1

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, mxArray const *prhs[])
{
    if (nrhs < NINPUTS - 1 || nlhs != NOUTPUTS)
        mexErrMsgTxt("must have at least NINPUTS - 1 input(s) and exactly one output");
    int maxT = (int) mxGetScalar(VAR_maxT);
    double np_d = mxGetScalar(VAR_np);
    int np = (int) np_d;
    int maxnsteps = (int) mxGetScalar(VAR_maxnsteps);
    int maxK = (int) mxGetScalar(VAR_maxK);
    int maxKp1substeps = (int) mxGetScalar(VAR_maxKp1substeps);
    int maxtotalsubsteps = (int) mxGetScalar(VAR_maxtotalsubsteps);
    int maxpresubsteps = (int) mxGetScalar(VAR_maxpresubsteps);
    if (((double) (MAX(maxpresubsteps, maxKp1substeps))) * np_d  > 2147483646.5)
        mexErrMsgTxt("Maximum allowed valued of max(maxpresubsteps, maxKp1substeps) * Nparticles is 2^31 - 1 = 2147483647");
    unsigned long long seedval;
    if (nrhs >= 8) {
        if (mxUINT64_CLASS != mxGetClassID(VAR_seedval) || mxGetNumberOfElements(VAR_seedval) != 1)
            mexErrMsgTxt("seedval must be a unit64 scalar");
        seedval = *((unsigned long long *) mxGetData(VAR_seedval));
    } else
        seedval = 1006;
    VAR_gdata = mxCreateNumericMatrix(1, sizeof(gpu5s_problem), mxUINT8_CLASS, mxREAL);
    gpu5s_problem* g = (gpu5s_problem *) mxGetData(VAR_gdata);
    g->h.fobs = NULL;
    g->h.u = NULL;
    int ecode = allocate_gpupfdata(g, maxT, np, maxnsteps, maxK, maxKp1substeps, maxtotalsubsteps, maxpresubsteps, seedval);
    
    if (ecode < 0) {
        mexPrintf("error code: %d\n", ecode);        
        mexErrMsgTxt("GPU memory allocation failed: error code: ");    
    }
}
