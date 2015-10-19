/*
Copyright (c) 2015, Nils Hammerla. n.hammerla@gmail.com
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include <math.h>
#include "mex.h"
#include "matrix.h"
#define MAX(X, Y) (((X) < (Y)) ? (Y) : (X))

/* 
Function to calculate sample entropy
Input:
  Data            double[]  Input sequence (double). Nx1 in matlab
  numSamples      int       Length of input
  wlen            int       Length of window (normally called m)
  r               double    Tolerance for "similarity"
  shift           int       Shift between samples (for subsampling). Shift=1 
                            corresponds to _no_ subsampling.

Output:
  double  -log (A / B), where A is the number of |D_m(i,j) < r|, and B is 
          |D_m+1(i,j) < r|. See wikipedia:
          https://en.wikipedia.org/wiki/Sample_entropy
*/
double calcSampleEntropy(double data[], int numSamples, int wlen, double r, int shift) {
    int A=0; 
    int B=0;
    int i,j,k;
    double m;
    
    /* ok, now we go through all windows data_i ... data_i+wlen and calculate the 
       Chebyshev distance to all the _following_ windows. As the distance is symmetric,
       d(i,j) = d(j,i), we only need to calculate the distance for less than half the pairs:
      
        D(i,j) (just calculate x)
            - x x x x
            - - x x x
            - - - x x
            - - - - x 
            - - - - -
    */ 
    for (i=0; i < numSamples-wlen*shift-shift; i += shift) {
        /* compare to all following windows > i */
        for (j=i+shift; j < numSamples-wlen*shift-shift; j+=shift) {
            m = 0; /* maximum so far */
            for (k=0; k < wlen; k++) 
                /* get max cheb. distance */
                m = MAX(m, fabs(data[i+k*shift]-data[j+k*shift]));
            /* first case, distance lower in first wlen positions */
            if (m < r) B++; 
            /* Second case, distance lower if we add the next element */
            if (MAX(m,fabs(data[i+wlen*shift]-data[j+wlen*shift])) < r) A++;
        }
    }
    /* return -log A/B */
    if (A>0 && B >0)
        return (-1 * log(((double) A) / ((double) B)));
    else
        return 0;
}

/* mex function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[]) {
    
    double result;
    result = calcSampleEntropy(mxGetPr(prhs[0]), (int) mxGetM(prhs[0]), (int) mxGetScalar(prhs[1]), (double) mxGetScalar(prhs[2]), (int) mxGetScalar(prhs[3]));
    plhs[0] = mxCreateDoubleScalar(result);
}


