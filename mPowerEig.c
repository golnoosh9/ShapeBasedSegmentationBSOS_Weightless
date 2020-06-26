#include <stdio.h>
#include <math.h>
#include "mex.h"

#define MAX 10

/* 
 * This is the power iteration method to find the maximum
 * eigenvalue/eigenvector a n-by-n matrix. This method doesn't
 * require the matrix to be Hermitian for the maximum eigenvalue/eigenvecor.
 * But it DOES require the matrix to be Hermitian for the minimum
 * eigenvalue/vector. This approximation method may be improved by
 * setting a tolerance (currently the iteration is controlled by the number
 * of iterations, MAX).
 * 
 * Example: c = [1 0.5 0.2;0.5 1 0.5; 0.2 0.5 1];
 * then [u,v] = mPowerEig(c,0) is to find the largest eigenvalue/vector 
 * and  [u,v] = mPowerEig(c,1) is to find the minimum eigenvalue/vector
 * 
 * Reference: G.H. Golub, C.F. Van Load, "Matrix Computation"
 */
 
void power(double **a, int n, double *z, double *q, double *lamsub)
{
    int i,j,k;
    double norm, lambda;
    
    for (i=0; i<MAX; i++){
       for (j=0;j<n;j++) {
            z[j] = 0.0;
            for (k=0; k<n;k++) z[j] += a[j][k]*q[k];
        }
      
        /* find the norm of vector z */
        norm = 0.0;
        for (j=0;j<n;j++) {
            norm = norm + z[j]*z[j];
        }
        norm = sqrt(norm);

        /* scale vector z by norm */
        for (j=0;j<n;j++) q[j] = z[j]/norm;
        lambda = 0.0;
        for (j=0;j<n;j++) lambda = lambda + q[j] * z[j];
    }
    *lamsub = lambda;
}

/* entrance routine from Matlab mex-function */
void powermethod(double *x, int n, int mode, double *y, double *w)
{
    double **a, **b,*q, lambda,lambdamin,*z;
    int i,j;
    
    /* convert Matlab's matrix to C-convention matrix stroage */
    a = (double **) malloc(n * sizeof(double *));
	if (a == NULL) {
        printf("matrix a does not have enogh memory\n");
        return;
    }
    a[0] = (double *) malloc(n  * n * sizeof(double));
    if (a[0] == NULL) {
        printf("matrix a does not have enogh memory\n");
        return;
    }
	for (i=1;i<n;i++) a[i] = a[i-1]+n;

	/* fill in C-stroage matrix elements */
	for (i=0;i<n;i++) {
	    for(j=0;j<n;j++) {
	        a[i][j] = x[i+n*j];
	    }
	    printf("\n");
	}    
    
    /* initialize the iterated eigenvector and working vector*/
    z = (double *) malloc(n*sizeof(double));
    q = (double *) malloc(n*sizeof(double));
    for (i=0;i<n;i++) q[i] = 0.0;
    q[0] = 1.0;
    
    /* find the maximum eigenvale/vector of matrix a */
    power(a,n,z,q,&lambda);
    
    /* mode ==0 is for maximum eigenvalue/vector */
    if (mode == 0) {
        for (i=0;i<n;i++) {
            y[i] = q[i];
        }
       *w = lambda;
       free(a);
       free(q);
       free(z);
    }
    /* else mode ~= 0, for minimum eigenvalue/vector */
    else {
        b = (double **) malloc(n * sizeof(double *));
	    if (b == NULL) {
           printf("matrix a does not have enogh memory\n");
           return;
        }
        b[0] = (double *) malloc(n  * n * sizeof(double));
        if (b[0] == NULL) {
            printf("matrix b does not have enogh memory\n");
            return;
        }    
    	for (i=1;i<n;i++) b[i] = b[i-1]+n;
     
    for (i=0;i<n;i++) {
        for(j=0;j<n;j++) b[i][j] = -a[i][j];
    }
    for (i=0;i<n;i++) b[i][i]=lambda-a[i][i];
    
    /* perform one more time power iteration method */
    	for (i=0;i<n;i++) q[i] = 0.0;
        q[0] = 1.0;
            
    /* iterative method to find maximum eigenvalue */
    power(b,n,z,q,&lambdamin);
    
    for (i=0;i<n;i++) y[i]=-q[i];
    *w = lambda - lambdamin;
    
       free(a);
       free(q);
       free(z);
       free(b);
    } 
}

/* the gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
  int n, mode;
  double *x, *y,*z;
  
  mexPrintf("Usage: [u,v] = mPowerEig(a,mode)\n");
  mexPrintf("where: a is the input matrix, n-by-n\n");
  mexPrintf("           mode = 0 is for the maximum eigenvalue/eigenvector\n");
  mexPrintf("                     =1 is for the minimum eigenvalue/eigenvector\n");
  mexPrintf("            u is the max/min eigenvector\n");
  mexPrintf("            v is the max/min eigenvalue\n");
  
  /* get the dimension of the data */
  x = mxGetPr(prhs[0]);
  n = mxGetN(prhs[0]);
  mode = mxGetScalar(prhs[1]);
  
  /*  set the output pointer to the output matrix */
  plhs[0] = mxCreateDoubleMatrix(n,1,mxREAL);
  plhs[1] = mxCreateDoubleScalar(0.0);
  y = mxGetPr(plhs[0]);
  z = mxGetPr(plhs[1]);
  
  /*  call the C subroutine */
  powermethod(x,n,(int)mode,y,z);
}
