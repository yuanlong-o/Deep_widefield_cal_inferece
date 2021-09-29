/* path_out = dendrite_randomwalk_cpp(float *M, int *root, int *ends, 
 *   float distsc, int maxlength, float fillweight, int maxel,  
 *   int minlength, int *volsize, int *path, int *bglengthp)
 *
 * Dendrite random walk for cortical dendrite simulation, mex source
 *
 *
 * 
 *  - M          - matrix that indicates the difficulty to occupy a
 *                 particular location
 *  - root       - 1x3 start location for the random walk path
 *  - ends       - Nx3 directed ending locations for the random walk
 *  - distsc     - weighting parameter greater than 0 setting relative
 *                 weight of directed walk
 *  - maxlength  - maximum length for the random path
 *  - fillweight - weighting value for a single step against occupancy rate
 *  - maxel      - maximum number of components within a single voxel
 *  - minlength  - minimum length for the random path
 *
 *  - path_out   - ourput paths of the random walks
 * 
 * 2017 Alex Song (as31@princeton.edu)
 */

#include <mex.h>
#include <math.h>
#include <vector>
using namespace std;

void dendrite_randomwalk_cpp(float *M, int *root, int *ends, float distsc, int maxlength, float fillweight, int maxel, int minlength, int *volsize, int *path, int *bglengthp)
{
  int i,j;
  int jmin;
  int *curr = root;
  int currIdx,endsIdx;
  vector<float> Mvals (maxlength,0);
  float dist,minmat,maxfill;
  vector<float> matvals (6,0);
  vector<float> distvec (3,0);
  currIdx = curr[0]+curr[1]*volsize[0]+curr[2]*volsize[0]*volsize[1];
  endsIdx = ends[0]+ends[1]*volsize[0]+ends[2]*volsize[0]*volsize[1];
  int bglength = 0;
  maxfill = maxel*fillweight;
  for(i = 0; i<maxlength; i++){
    dist = sqrt(pow((float)(ends[0]-curr[0]),2)+pow((float)(ends[1]-curr[1]),2)+pow((float)(ends[2]-curr[2]),2))/distsc;
    distvec[0] = (curr[0]-ends[0])/dist;
    distvec[1] = (curr[1]-ends[1])/dist;
    distvec[2] = (curr[2]-ends[2])/dist;
    
    jmin = 6;
    minmat = FLT_MAX;
    
    if(curr[0]!=(volsize[0]-1)){
      matvals[0] = M[currIdx+1]+distvec[0];
      if(matvals[1]<minmat){
        jmin = 0;
        minmat = matvals[0];
      }
    }
    if(curr[1]!=(volsize[1]-1)){
      matvals[1] = M[currIdx+volsize[0]]+distvec[1];
      if(matvals[1]<minmat){
        jmin = 1;
        minmat = matvals[1];
      }
    }
    if(curr[2]!=(volsize[2]-1)){
      matvals[2] = M[currIdx+volsize[0]*volsize[1]]+distvec[2];
      if(matvals[2]<minmat){
        jmin = 2;
        minmat = matvals[2];
      }
    }
    if(curr[0]!=0){
      matvals[3] = M[currIdx-1]-distvec[0];
      if(matvals[3]<minmat){
        jmin = 3;
        minmat = matvals[3];
      }
    }
    if(curr[1]!=0){
      matvals[4] = M[currIdx-volsize[0]]-distvec[1];
      if(matvals[4]<minmat){
        jmin = 4;
        minmat = matvals[4];
      }
    }
    if(curr[2]!=0){
      matvals[5] = M[currIdx-volsize[0]*volsize[1]]-distvec[2];
      if(matvals[5]<minmat){
        jmin = 5;
        minmat = matvals[5];
      }
    }
    
    if(minmat<maxfill){
      switch(jmin) {
        case 0: curr[0]++; break;
        case 1: curr[1]++; break;
        case 2: curr[2]++; break;
        case 3: curr[0]--; break;
        case 4: curr[1]--; break;
        case 5: curr[2]--; break;
      }
      currIdx = curr[0]+curr[1]*volsize[0]+curr[2]*volsize[0]*volsize[1];
      path[i*3] = curr[0];
      path[i*3+1] = curr[1];
      path[i*3+2] = curr[2];
      Mvals[i] = minmat;
      M[currIdx] = FLT_MAX;
      if(curr[0] == 0 || curr[1] == 0 || curr[2] == 0 || curr[0] == (volsize[0]-1) || curr[1] == (volsize[1]-1) || curr[2] == (volsize[2]-1)){
        break;
      }
    } else {
      bglength = i-1;
      break;
    }
    if(currIdx == endsIdx){
      break;
    }
  }
  if(bglength <= 0){
    bglength = i;
  }
  if(bglength >= minlength){
    for(i = 0; i<bglength; i++){
      if(Mvals[i] < maxfill) {
        M[path[i*3]+path[i*3+1]*volsize[0]+path[i*3+2]*volsize[0]*volsize[1]] = Mvals[i]+fillweight;
      }
    }
  } else {
    for(i = 0; i<bglength; i++){
      M[path[i*3]+path[i*3+1]*volsize[0]+path[i*3+2]*volsize[0]*volsize[1]] = Mvals[i];
    }    
  }
  bglengthp[0] = bglength;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  float *M;
  int *root;
  int *ends;
  float distsc;
  int maxlength, minlength;
  float fillweight;
  int maxel;
  int *path;
  int *pathc;
  int *bglength;
  int root2[3];
  int ends2[3];
  int i;
  
  int *volsize;  
  volsize = (int*)mxGetData(prhs[8]);
            
  M = (float*)mxGetData(prhs[0]);
  root = (int*)mxGetData(prhs[1]);
  ends = (int*)mxGetData(prhs[2]);

  root2[0] = root[0]-1;
  root2[1] = root[1]-1;
  root2[2] = root[2]-1;
  ends2[0] = ends[0]-1;
  ends2[1] = ends[1]-1;
  ends2[2] = ends[2]-1;
  
  distsc = (float)mxGetScalar(prhs[3]);
  maxlength = (int)mxGetScalar(prhs[4]);
  fillweight = (float)mxGetScalar(prhs[5]);
  maxel = (int)mxGetScalar(prhs[6]);
  minlength = (int)mxGetScalar(prhs[7]);

  plhs[0] = mxCreateNumericMatrix(3, maxlength, mxINT32_CLASS, mxREAL);
  path = (int*)mxGetData(plhs[0]);
  plhs[1] = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
  bglength = (int*)mxGetData(plhs[1]);

  dendrite_randomwalk_cpp(M,root2,ends2,distsc,maxlength,fillweight,maxel,minlength,volsize,path,bglength);

  plhs[2] = mxCreateNumericMatrix(*bglength, 3, mxINT32_CLASS, mxREAL);
  pathc = (int*)mxGetData(plhs[2]);
  
  for(i = 0; i<(*bglength); i++){
    pathc[i] = path[i*3]+1;
    pathc[i+(*bglength)] = path[i*3+1]+1;
    pathc[i+2*(*bglength)] = path[i*3+2]+1;
  }

}