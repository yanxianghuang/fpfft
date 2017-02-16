//#########################################################
//# Project         : Fixed-point FFT
//#
//#
//# Author          : Yanxiang Huang
//#
//# Company         : imec
//#                   Kapeldreef 75
//#                   B-3001 Leuven
//#                   Belgium
//#                   http://www.imec.be
//##########################################################

// input format:
// arg1-4: real_in, imag_in, real_tw, imag_tw
// arg5-8: Number_of_FFT_points, Number_of_FFT_stages, Input_frac_size, tw_frac_size;

#include "mex.h"
//#include <string.h> //for memcpy/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <cstdint.h> // use this for pure cpp
#include <tmwtypes.h>

//#define DEBUG


// decimation-in-time FFT
void my_fix_fft(int64_T *x_real, int64_T *x_imag, int32_T *tw_real, int32_T *tw_imag, int N, int Nstage, int tw_frac ){// logN is base 2 log(N){
  
  int innerSize = N;
  int Nblocks   = 1;
  int in1, in2, bf;
  int stage, block_space, block, inner;
  
  tw_frac = tw_frac+1; // simplified computation, see [shift-right-tw]
  
  int64_T temp_real, temp_imag;
  
  for(stage = 0; stage < Nstage; stage++){// stage
    
    block_space = innerSize;
    innerSize = (innerSize>>1);
    
    #ifdef DEBUG
    printf("DEBUG, inner size changde to %d\n", innerSize);
    #endif
    
    for(block = 0; block < Nblocks; block++){ 
      
      #ifdef DEBUG
      printf("DEBUG, block is now %d\n", block);
      #endif
      
      for(inner = 0; inner < innerSize; inner++){
        
        // index for butterfly, and for bf coeffs
        in1 = (block*block_space) + inner;
        in2 = in1 + innerSize;
        bf  = inner*Nblocks;
        
        // butter-fly
        temp_real = x_real[in1] - x_real[in2];
        temp_imag = x_imag[in1] - x_imag[in2];       
        x_real[in1] = (x_real[in1] + x_real[in2]) >> 1; // Note that for additions, right-shift by one, so the wordlength of data is always constant
        x_imag[in1] = (x_imag[in1] + x_imag[in2]) >> 1;      
        x_real[in2] = (temp_real*tw_real[bf] - temp_imag*tw_imag[bf]) >> tw_frac; // Note [shift-right-tw],  shift-right-by the frac-width of tw, keeping the wordlength constant.
        x_imag[in2] = (temp_real*tw_imag[bf] + temp_imag*tw_real[bf]) >> tw_frac;
        
        #ifdef DEBUG
        printf("%d, %d, %d\n", in1, in2, bf);
        printf("value: %lf, %lf*i;   %lf, %lf*i\n", ((double)x_real[in1])/(1<<30), ((double)x_imag[in1])/(1<<30),
          ((double)x_real[in2])/(1<<30), ((double)x_imag[in2])/(1<<30));
        #endif
      }// end of block 
    }// end iteration
    Nblocks  = (Nblocks <<1);
  }// end FFT
  return;
}







void mexFunction( int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
  //checks
  if (nrhs != 4+4){
    mexErrMsgTxt("4 double inputs required: real (double); imag (double), real_tw (double), imag_tw (double),\n 6 int option required: fft_size,  fft_stage, x_frac, tw_frac");
  }
  if ( !mxIsDouble(prhs[0]) || !mxIsDouble(prhs[1]) || !mxIsDouble(prhs[2]) || !mxIsDouble(prhs[3]) ){
    mexErrMsgTxt("first 4 input should be double");
  }
  
  
  // double inputs
  double *real;
  double *imag;
  double *real_tw;
  double *imag_tw;
  real = mxGetPr(prhs[0]);
  imag = mxGetPr(prhs[1]);
  real_tw = mxGetPr(prhs[2]);
  imag_tw = mxGetPr(prhs[3]);
  
  // int inputs
  int fft_size, fft_stage, x_frac, tw_frac;
  fft_size  =  (int ) mxGetScalar(prhs[4]);
  fft_stage =  (int ) mxGetScalar(prhs[5]);
  x_frac    =  (int ) mxGetScalar(prhs[6]);
  tw_frac   =  (int ) mxGetScalar(prhs[7]);
  
  
  // cast it to fixed-point numbers
  // fixed-size
  int64_T fi_real [2048];
  int64_T fi_imag [2048];
  int32_T fi_real_tw [2048];
  int32_T fi_imag_tw [2048];
  //int64_T *fi_real = new int64_T[fft_size];
  //int64_T *fi_imag = new int64_T[fft_size];
  //int32_T *fi_real_tw = new int32_T[fft_size/2];
  //int32_T *fi_imag_tw = new int32_T[fft_size/2];
  int i;
  for(i=0;i<fft_size;i++){
    fi_real[i] = (int64_T) ( real[i] *((1<<x_frac)) );
    fi_imag[i] = (int64_T) ( imag[i] *((1<<x_frac)) );
    //printf("%lld; \n", fi_real[i]);
  }
  
  for(i=0;i<fft_size/2;i++){
    fi_real_tw[i] = (int32_T) (real_tw[i]*(1<<tw_frac));
    fi_imag_tw[i] = (int32_T) (imag_tw[i]*(1<<tw_frac));
    //printf("%ld; \n", fi_real_tw[i]);
  }
  
  //printf("%d, %d, %d, %d\n", sizeof(short), sizeof(int), sizeof(long), sizeof(long long));
  //printf("%d, %d\n", sizeof(int32_T), sizeof(int64_T));
  
  
  // FFT
  my_fix_fft(fi_real, fi_imag, fi_real_tw, fi_imag_tw, fft_size, fft_stage, tw_frac);
  
  //output
  double *dp1, *dp2;
  plhs[0] = mxCreateDoubleMatrix(1,fft_size,mxREAL);
  plhs[1] = mxCreateDoubleMatrix(1,fft_size,mxREAL);
  
  dp1 = mxGetPr(plhs[0]);
  dp2 = mxGetPr(plhs[1]);
  // cast to double
  for(i=0; i<fft_size; i++){    
    dp1[i] =   ((double) fi_real[i]) * ((double)(1<<fft_stage)) / ((double)(1<<x_frac));
    dp2[i] =   ((double) fi_imag[i]) * ((double)(1<<fft_stage)) / ((double)(1<<x_frac));
  }
  //printf("new\n");
  
  //free-up
  //free(fi_real);free(fi_imag); free(fi_real_tw); free(fi_imag_tw);
}
