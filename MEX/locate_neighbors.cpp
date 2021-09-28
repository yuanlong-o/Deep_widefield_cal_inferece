#include <math.h>
#include "mex.h"
#include "matrix.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (!mxIsSingle(prhs[0])){    
        mexErrMsgTxt("Array must be single.");                             // Make sure first input is a single
    }
    if (!mxIsSingle(prhs[1])){
        mexErrMsgTxt("Array must be single.");                             // Make sure second input is a single
    }
    if (!mxIsSingle(prhs[2])){
        mexErrMsgTxt("Array must be single.");                             // Make sure third input is a single
    }
    if (!mxIsSingle(prhs[3])){    
        mexErrMsgTxt("Array must be single.");                             // Make sure fourth input is a single
    }
    if (!mxIsSingle(prhs[4])){
        mexErrMsgTxt("Array must be single.");                             // Make sure fifth input is a single
    }
    float *X = (float*)mxGetPr(prhs[0]);                                   // Declare pointer to X-axis locations
    float *Y = (float*)mxGetPr(prhs[1]);                                   // Declare pointer to Y-axis locations
    float *Z = (float*)mxGetPr(prhs[2]);                                   // Declare pointer to Z-axis locations
    float *ctr = (float*)mxGetPr(prhs[3]);                                 // Declare pointer to center location
    float *thresh = (float*)mxGetPr(prhs[4]);                              // Declare pointer to the threshold
    long int ii = 0;
    size_t Nx = mxGetNumberOfElements(prhs[0]);                            // Get number of elements for the X-locations
    size_t Ny = mxGetNumberOfElements(prhs[1]);                            // Get number of elements for the Y-locations
    size_t Nz = mxGetNumberOfElements(prhs[2]);                            // Get number of elements for the Z-locations
    if ( Nx!=Ny || Nx!=Ny || Ny!=Nz){
        mexErrMsgTxt("Number of elements in X, Y, and Z must mach.");      // Make sure number of elements match
    }
    size_t Npt = mxGetNumberOfElements(prhs[3]);                           // Get number of elements for the center location
    if (Npt != 3){
        mexErrMsgTxt("Need a 3-vector as the center point.");              // Check that the location is a 3-vector
    }
    plhs[0] = mxCreateLogicalArray(3, mxGetDimensions(prhs[0]));           // Create logical output array of same size as X-, Y-, and Z-location arrays
    mxLogical *V = mxGetLogicals(plhs[0]);                                 // Make a pointer to the output
    for (ii=0; ii<Nx; ii++) {        
        V[ii] = (mxLogical*)(sqrt(pow(X[ii] - ctr[0],2) + pow(Y[ii] - ctr[1],2) + pow(Z[ii] - ctr[2],2)) < thresh[0]); // Calculate if the iith point is within thresh from the center location
    }
}


