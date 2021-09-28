/*
 * Dendrite dijkstra for  cortical dendrite simulation, mex source
 * 
 * C++ code to grow dendrites via a stochastic process. This code
 * is the mex sourse that should be compiked through the latest
 * version of gpp through MATLAB in order to use with TAO-SIM. 
 * 
 * Details of the algorithm can be found in **
 * 
 * Alex Song (as31@princeton.edu)
 */

#include <mex.h>
#include <queue>
#include <vector>
using namespace std;

class compDist
{
public:
    bool operator()(pair<float,int> n1,pair<float,int> n2) {
        return n1.first>n2.first;
    }
};

// structure for queue
//typedef pair <int, float> edge;
//struct edge { int vert; float dist; };

void dendrite_dijkstra_cpp(float *M, int *pe, int root, int pdims, float *distance, int *pathfrom)
{
  int i;
  vector<bool> tovisit (pdims,true);
  priority_queue< pair<float, int>, vector<pair<float,int>>, compDist > unvisited;
  unvisited.push( make_pair(0,root) );
  fill_n(distance, pdims, FLT_MAX);
  distance[root] = 0;
  fill_n(pathfrom, pdims, INT_MAX);
  pair<float, int> cn;
  int nn;
  float ndist;
  while(!unvisited.empty()){
    cn = unvisited.top();
    unvisited.pop();
    for (i = 0; i< 6; i++){
      nn = cn.second+pe[i];
      if(nn>0 && nn<pdims){
        if(tovisit[nn]){
          ndist = cn.first+M[nn+pdims*i];
          if(ndist < distance[nn]){
            unvisited.push( make_pair(ndist,nn) );
            distance[nn] = ndist;
            pathfrom[nn] = cn.second;
          }
        }
      }
    }
    tovisit[cn.second] = 0;
  }
  for(i = 0; i<pdims; i++){
    pathfrom[i] += 1;
  }
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  float *M;
  int *pe;
  int root;
  int pdims;
          
  float *distance;
  int *pathfrom;
          
  M = (float*)mxGetData(prhs[0]);
  pe = (int*)mxGetData(prhs[1]);
  root = mxGetScalar(prhs[2])-1;
  pdims = mxGetM(prhs[0]);
  
  plhs[0] = mxCreateNumericMatrix(1, pdims, mxSINGLE_CLASS, mxREAL);
  distance = (float*)mxGetData(plhs[0]);

  plhs[1] = mxCreateNumericMatrix(1, pdims, mxINT32_CLASS, mxREAL);
  pathfrom = (int*)mxGetData(plhs[1]);

  dendrite_dijkstra_cpp(M,pe,root,pdims,distance,pathfrom);
}

// EOF
