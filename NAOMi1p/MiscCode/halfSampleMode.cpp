/**

http://www.sciencedirect.com/science/article/pii/S0167947305001581

*/


#include <mex.h>
#include <vector>
#include <algorithm>
#include <cmath>

<<<<<<< HEAD
#include "classes/MonotCubicInterpolator.hpp"
#include "classes/MonotCubicInterpolator.cpp"
=======
#include "MonotCubicInterpolator.hpp"
#include "MonotCubicInterpolator.cpp"
>>>>>>> 405214c411980eac8d04e3b2baa9f89c37090a29



//=============================================================================
template<typename Number>
Number findSampleMode(size_t n, const Number* x)
{
  while (true)
  {
    switch (n) {
    case 0:
      mexErrMsgIdAndTxt("halfSampleMode:sanity", "Should be impossible.");

    case 1:
      return x[0];

    case 2:
      return (x[0] + x[1])/2;

    case 3:
    {
      // Safety against unsigned data
      const double      diff              = ( double(x[1]) - double(x[0]) )
        - ( double(x[2]) - double(x[1]) )
        ;
      if (diff < 0)
        return (x[0] + x[1]) / 2;
      else if (diff > 0)
        return (x[1] + x[2]) / 2;
      else
        return x[1];
    }

    default:
    {
      const size_t      N                 = std::ceil( 0.5 * n );
      double            wMin              = std::numeric_limits<double>::infinity();
      size_t            j                 = n;
      for (size_t i = 0; i < n - N; ++i) {
        const double    width             = double(x[i + N - 1]) - double(x[i]);
        if (width < wMin) {
          wMin          = width;
          j             = i;
        }
      }
      x                += j;
      n                 = N;
    }
    }
  }
}


std::vector<size_t> chunkIndices(size_t numTotal, size_t numPerChunk)
{
  size_t                numChunks     = std::max<size_t>(1u, static_cast<size_t>(std::round(numTotal / numPerChunk)));
  std::vector<size_t>   chunk(numChunks + 1);

  const size_t          delta         = static_cast<size_t>(std::round(1.0 * numTotal / numChunks));
  for (size_t iChunk = 1; iChunk < numChunks; ++iChunk)
    chunk[iChunk]       = iChunk * delta;
  chunk[0]              = 0;
  chunk[numChunks]      = numTotal;

  return chunk;
}

template<typename Number, typename Ptr>
void sortFiniteItems(std::vector<Number>& sorted, Ptr first, Ptr last)
{
  sorted.clear();
  for (Ptr index = first; index != last; ++index)
    if (mxIsFinite(*index))
      sorted.push_back(*index);

  std::sort(sorted.begin(), sorted.end());
}


//=============================================================================
template<typename Number>
void findSampleMode(const mxArray* matData, mxArray*& mode, size_t window)
{
  const Number*         data              = (const Number*) mxGetData(matData);
  const size_t          nData             = mxGetM(matData);
  const size_t          nPixels           = mxGetN(matData);
  if (nData < 1)
    mexErrMsgIdAndTxt("halfSampleMode:data", "data must contain at least one item.");


  mode                  = mxCreateNumericMatrix(window ? nData : 1, nPixels, mxGetClassID(matData), mxREAL);
  Number*               modePtr           = (Number*) mxGetData(mode);
  std::vector<Number>   sortedData;
  sortedData.reserve(nData);

  if (window) {

    std::vector<size_t> chunks            = chunkIndices(nData, window);
    const size_t        numChunks         = chunks.size() - 1;

    std::vector<double> yChunk(numChunks);
    std::vector<double> xNominal(numChunks), xTemp(numChunks);
    for (size_t iChunk = 0; iChunk < numChunks; ++iChunk)
      xNominal[iChunk]  = (chunks[iChunk] + chunks[iChunk + 1]) / 2.0;

      
    Number*             pixMode           = modePtr;
    for (size_t iPix = 0; iPix < nPixels; ++iPix, data += nData) {
      // Compute mode in chunks
      bool              missingData       = false;
      for (size_t iChunk = 0; iChunk < numChunks; ++iChunk) {
        sortFiniteItems(sortedData, data + chunks[iChunk], data + chunks[iChunk+1]);
        if (sortedData.empty()) {
          yChunk[iChunk]= mxGetNaN();
          missingData   = true;
        }
        else            yChunk[iChunk]    = findSampleMode(sortedData.size(), sortedData.data());
      }

      // Handle cases of incomplete data
      std::vector<double>*    xChunk;
      if (missingData) {
        xTemp.clear();
        for (size_t iChunk = 0; iChunk < numChunks; ++iChunk) {
          if (!mxIsNaN(yChunk[iChunk])) {
            yChunk[xTemp.size()]  = yChunk[iChunk];
            xTemp.push_back(xNominal[iChunk]);
          }
        }
        xChunk          = &xTemp;
        yChunk.resize(xChunk->size());
      } 
      else xChunk       = &xNominal;


      // Interpolate between chunks
      Opm::MonotCubicInterpolator         interp(*xChunk, yChunk);
      for (size_t iData = 0; iData < nData; ++iData, ++pixMode)
        *pixMode        = interp.evaluate(iData);

      if (missingData)
        yChunk.resize(numChunks);
    }

  }

  else {
    for (size_t iPix = 0; iPix < nPixels; ++iPix, data += nData) {
      sortFiniteItems(sortedData, data, data + nData);
      modePtr[iPix]     = findSampleMode(sortedData.size(), sortedData.data());
    }
  }
}


//=============================================================================
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  //----- Parse arguments
  if (nrhs < 1 || nrhs > 2 || nlhs > 1) {
    mexErrMsgIdAndTxt ( "halfSampleMode:usage"
                      , "Usage:     mode = halfSampleMode(data, [window = inf])\n"
                        "All matrices are expected to n by p where n is the number of samples per column dataset."
                      );
  }

  const mxArray*        data            = prhs[0];
  size_t                window          = 0;
  if (nrhs > 1 && !mxIsEmpty(prhs[1])) {
    const double        input           = mxGetScalar(prhs[1]);
    if (mxIsFinite(input) && input < mxGetM(data))
      window            = static_cast<size_t>(input);
    else
      window            = mxGetM(data);
  }


  switch (mxGetClassID(data)) {
  case mxSINGLE_CLASS :   findSampleMode<float         >(data, plhs[0], window);   break;
  case mxCHAR_CLASS   :   findSampleMode<char          >(data, plhs[0], window);   break;
  case mxDOUBLE_CLASS :   findSampleMode<double        >(data, plhs[0], window);   break;
  case mxINT8_CLASS   :   findSampleMode<char          >(data, plhs[0], window);   break;
  case mxUINT8_CLASS  :   findSampleMode<unsigned char >(data, plhs[0], window);   break;
  case mxINT16_CLASS  :   findSampleMode<short         >(data, plhs[0], window);   break;
  case mxUINT16_CLASS :   findSampleMode<unsigned short>(data, plhs[0], window);   break;
  case mxINT32_CLASS  :   findSampleMode<int           >(data, plhs[0], window);   break;
  case mxUINT32_CLASS :   findSampleMode<unsigned int  >(data, plhs[0], window);   break;
  case mxINT64_CLASS  :   findSampleMode<int64_t       >(data, plhs[0], window);   break;
  case mxUINT64_CLASS :   findSampleMode<uint64_t      >(data, plhs[0], window);   break;

  default:
    mexErrMsgIdAndTxt("halfSampleMode:arguments", "Unsupported type of data.");
  }
}

