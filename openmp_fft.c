//
//  serial_fft.c
//  
//
//  Created by MARELAS GEORGE on 1/15/15.
//
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>
#include <time.h>
#include <omp.h>

#define N 11 //power of for samples 
#define numthreads 16//choose # of threads

void fft(double *x, double *y, int n);
int reverseBits(int num);

int main(int argc, char **argv) {
    
    //int N =2;
    int n = 0;
    int i = 0;
    double *x;
    double *y;
    clock_t t;
    double seconds;
    
    
    srand(time(NULL));
    n = pow(2,N);
    printf("Samples are 2^ %d = %d",N,n);
    
    x = (double*)calloc(n,sizeof(double));
    if(x==NULL){
        printf("no mem allocate for x samples\n");
    }
    
    y = (double*)calloc(n,sizeof(double));
    if(y==NULL){
        printf("no mem allocate for y outuput\n");
            
    }
    
    
    for( i=0; i<n;  i++  ) {
        x[i] = ((rand()%3)+1)/0.3;
    }
    
    t = clock(); 
    fft(x, y, n);
    t = clock() - t;
    seconds = ((double)t)/CLOCKS_PER_SEC;
    printf("\n%f sec for n = %d and threads = %d\n",seconds,n,numthreads);
    free(x);
    free(y);
    
    
    return ( 0 );
}


void fft(double *x, double *y, int n) {
    
    double complex w;
    
    int i,j,m;
    double y_elements = 0;
    int temp_j = 0xffffe;
    int temp_k = 0x00000001;
    int k;
    int ii = 0;
    FILE *fp;
    fp = fopen("results.txt","w");
    
    int tid,nthreads,chunk;
    
    
    double *R = NULL;
    double *S  =NULL;
    
    w = cexp((-2*M_PI*I)/n);
    printf("The w is: %lf + ( %lfi )\n", creal(w), cimag(w));
    
    y_elements  = pow(creal(w),i);
    printf("%.2lf\n", y_elements);
    
    
    R = (double*)calloc(n,sizeof(double));
    if(R==NULL){
        printf("no mem allocate for R\n");
        return;
    }
    S = (double*)calloc(n,sizeof(double));
    if(S==NULL){
        printf("no mem allocate for S\n");
        return;
    }
    
    for (i = 0; i<n; i++) {
        R[i] = x[i];
    }
    
    #if 0
    for (i = 0; i<n; i++) {
      fprintf(fp,"X[%d] =  %.2lf => R[%d] = %.2lf\n",i, x[i], i, R[i]);
    }
    #endif

    
      //#pragma omp parallel for //schedule(dynamic,chunk) //shared(tid) private(x,y,S,R,n)
      for (m = 0;  m<N;  m++) {
	  for (i = 0;  i<n; i++) {
	      S[i] = R[i];
	  }
        
	  i = 0;
	  //printf("tempk: %d\n", );
   
	  omp_set_num_threads(numthreads);
	  chunk = n/numthreads;
	  #pragma omp parallel private(tid,i,j,k) 
	  {
	      tid = omp_get_thread_num();
	      //printf("thread number = %d\n",tid);
	      // Only master thread does this 
	      if (tid == 0) 
	      {
		  nthreads = omp_get_num_threads();
		 // printf("Number of threads = %d\n", nthreads);
	      }
	      // #pragma omp parallel for schedule(dynamic,chunk)
	      for (i = tid*chunk; i<(tid+1)*chunk; i++) {
              
		j = i & temp_j;
		k = i | temp_k;
            
		R[i] = S[j] + S[k]*( pow(creal(w), i) );       
		}
	    }  
	  //printf("%d\n", temp_k );
        
	  temp_k = temp_k << 1;
	  temp_j = temp_j << 1;
	  temp_j = temp_j | 0x1;
	
    }
        


    for (i = 0; i<n; i++) {
      int NO_OF_BITS = sizeof(i)*8;
      //printf("**%d**",NO_OF_BITS);
      fprintf(fp,"R[%d] = %.2lf\n", i, R[i]);
      
    }
    


    for (i = 0; i<n; i++) {
        
	int reverse = reverseBits(i);  
	y[i] = R[reverse];
	fprintf(fp,"Y[%d] = %.2lf\n", i, y[i]);
    
    }
 
    
    free(R);
    free(S);   
    
}

int reverseBits(int num){
    {
    int  NO_OF_BITS = N;
    int reverse_num = 0, i, temp;
 
    for (i = 0; i < NO_OF_BITS; i++)
    {
        temp = (num&(1 << i));
        if(temp==1)
            reverse_num |= (1 << ((NO_OF_BITS - 1) - i));
    }
 
    return reverse_num;
}
  
}
