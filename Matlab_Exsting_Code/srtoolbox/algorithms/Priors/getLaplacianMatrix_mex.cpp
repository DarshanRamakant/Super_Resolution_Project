#include <mex.h>
#include <math.h>
#include <vector>

using namespace std;

typedef struct
{
    double x;
    double y;
} Pixel;

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
vector<Pixel> CreatePixelSequence(int dimY, int dimX);
mxArray* AllocateLoGMatrix(int n, mwSize* nzmax);
mxArray* ComposeLoGMatrix(const vector<Pixel>& pixelVec);

/**
    MEX implementation for the creation of the LoG matrix D such that D*x
    represents the LoG (Laplacian of Gaussian) of x.
    
    - prhs contains the image size in format (y, x) for an image in matrix
      notation where the LoG is computed.
    - plhs contains the sparse LoG matrix.
*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{ 
    // Get dimensions and determine 2-D coordinates of pixels that are
    // organized in a linear order (1,1), (1,2), (1,3), ... (m,n)
    int dimY = mxGetScalar(prhs[0]);
    int dimX = mxGetScalar(prhs[1]);
    vector<Pixel> pixelVec = CreatePixelSequence(dimY, dimX);
    
    plhs[0] = ComposeLoGMatrix(pixelVec);
}

/**
    Determine the 2-D coordinates of pixels from an image with specified
    dimensions from a linearized order of these pixels.
*/
vector<Pixel> CreatePixelSequence(int dimY, int dimX)
{
    vector<Pixel> pixelVec;
    for (int y = 0; y < dimY; y++)
    {
        for (int x = 0; x < dimX; x++)
        {
            Pixel p;
            p.x = x;
            p.y = y;
            pixelVec.push_back(p);
        }  
    }
    return pixelVec;
}

/**
    Compose the sparse LoG matrix.
*/
mxArray* ComposeLoGMatrix(const vector<Pixel>& pixelVec)
{
    mwSize nzmax;
    int n = pixelVec.size();
    
    // Allocate memory for the sparse matrix and get pointers to matrix data.
    mxArray* logMat = AllocateLoGMatrix(n, &nzmax);
    double* sr  = mxGetPr(logMat);
    mwIndex* irs = mxGetIr(logMat);
    mwIndex* jcs = mxGetJc(logMat);
   
    mwIndex k = 0;
    for (mwSize j = 0; j < n; j++)     // Iterate over all columns
    {
        jcs[j] = k;
        for (mwSize i = 0; i < n; i++)     // Iterate over all rows
        {
            double weight = 0;
            if ( ( abs(pixelVec[i].x - pixelVec[j].x) + abs(pixelVec[i].y - pixelVec[j].y) ) == 0 )
            {    
                weight = 1;
            }
            else if ( ( abs(pixelVec[i].x - pixelVec[j].x) + abs(pixelVec[i].y - pixelVec[j].y) ) == 1 )
            {
                weight = -0.25;
            }
            
            if (weight != 0)
            { 
                sr[k] = weight;
                irs[k] = i;     // We have found a row with an non-zero element.     
                k++;
            }
        }
  }
  jcs[n] = k;
  
  return logMat;
}

/**
    Allocate and initialize the sparse LoG matrix.
*/
mxArray* AllocateLoGMatrix(int n, mwSize* nzmax)
{
    // Allocate space for sparse matrix 
    double percent_sparse = (15.0*n) / (1.0*n*n);
    mwSize nz = (mwSize)ceil((double)n*(double)n*percent_sparse);
    if (nz < 10)
    {
        nz = 10;
    }
    mxArray* logMat = mxCreateSparse(n, n, nz, mxREAL);
    *nzmax = nz;
    return logMat;
}