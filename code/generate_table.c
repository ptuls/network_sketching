/* ============================================================= */
/* === MATLAB/generate_table mexFunction ======================= */
/* ============================================================= */
/* ----------------------------------------------------------------
 * T = generate_table(l,b)
 *
 * where:
 * INPUT
 *
 *    l: level
 *    b: number of bits in the sketch
 *
 * OUTPUT
 *
 *    T: generated table 
 *
 * Copyright (c) Paul Tune <paul.tune@adelaide.edu.au> 24 Jul 2014
 * ----------------------------------------------------------------
 */

#include <math.h>
#include <stdlib.h>  /* needed for qsort() */
#include "mex.h"
#include "matrix.h"

void mexFunction(
        int nlhs,       mxArray *plhs[],
        int nrhs, const mxArray *prhs[]
        )
{
  /* Declare variables */
  int i,j,k,l,L,G,b;
  int c,mask,e;
  double *lx,*bx,*Tx,**T;

  /* Check for proper number of input and output arguments. */    
  if (nrhs != 2) {
    mexErrMsgTxt("Two input arguments required: T = generate_table(l,b)");
  }
  if (nlhs > 1) {
    mexErrMsgTxt("Too many output arguments.");
  }

  /* Check data type of input argument. */
  /*if (mxGetNumberOfDimensions(prhs[0]) != 1 ||
      mxGetNumberOfDimensions(prhs[1]) != 1) {
    mexErrMsgTxt("First 2 input arguments must be one dimensional\n");
  }*/
  
  /* Get parameters */
  lx = mxGetPr(prhs[0]);
  bx = mxGetPr(prhs[1]);
  
  /* Other parameters */
  b = bx[0];
  l = lx[0];
  G = l*b;
  L = 1 << G;
  
  /* create output matrix */
  plhs[0] = mxCreateDoubleMatrix(L, L, mxREAL);
  Tx = mxGetPr(plhs[0]);

  T  = (double**) mxMalloc(L*sizeof(double*));
  for (i = 0, j = 0; i < L; i++, j+=L) {
    T[i] = Tx + j;
  }
  
  /* initialise array */
  for (i = 0; i < L*L; i++) {
    Tx[i] = 0.0;
  }
  
  mask = (1 << l)-1;
  for (i=0; i < L; i++){
      for (j=0; j < L; j++){
          e = i^j;
          c = parity(e & mask);
          for (k=1; k < b; k++){
             e = e >> l;
             c |= parity(e & mask) << k;
          }
          T[i][j] =  c;
      }
  }
                  
  /* garbage collection */
  mxFree(T);
}

int parity(v)
{
    v ^= v >> 16;
    v ^= v >> 8;
    v ^= v >> 4;
    v &= 0xf;
    return (0x6996 >> v) & 0x1;
}