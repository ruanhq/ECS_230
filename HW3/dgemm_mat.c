#include "stdio.h"
#include "stdlib.h"
#include "time.h"
#include "sys/time.h"
#include "math.h"
#include "cblas.h"
//References:
//1. https://stackoverflow.com/questions/39002052/how-i-can-print-to-stderr-in-c
//2. https://en.wikipedia.org/wiki/C_dynamic_memory_allocation
//3. https://stackoverflow.com/questions/28654438/matrix-vector-product-with-dgemm-dgemv
//4. https://www.tutorialspoint.com/c_standard_library/c_function_clock.htm
//5. https://canvas.ucdavis.edu/courses/246067/files/folder/Homework/timing1.c
//6. https://canvas.ucdavis.edu/courses/246067/files/folder/Homework/timing2.c
//7. http://www.netlib.org/blas


// BEGIN SUBparts
//Directly replicate Professor Gygi's code:gtod and readTSC 
volatile double gtod(void) {
  /* Directly from timing1.c on ~fgygi */
  static struct timeval tv;
  static struct timezone tz;
  gettimeofday(&tv,&tz);
  return tv.tv_sec + 1.e-6*tv.tv_usec;
}

long long readTSC(void)
{
  union { long long complete; unsigned int part[2]; } ticks;
  __asm__ ("rdtsc; mov %%eax,%0;mov %%edx,%1"
    : "=mr" (ticks.part[0]), "=mr" (ticks.part[1])
    : /* no inputs */
    : "eax", "edx");
  return ticks.complete;
}
// END SUBROUTINES

// BEGIN MAIN
int main(int argc, char** argv)
{
  if(argc<2){
    exit(1);
  }
  int n_row=atoi(argv[1]);
  int n_rep=atoi(argv[2]); 
  int i,j,k,l; 
  double sum=0.0;
  // dynamic memory allocation for the matrices
  double *matA=(double*)malloc((n_row*n_row)* sizeof(double));
  double *matB=(double*)malloc((n_row*n_row)* sizeof(double));
  double *matC=(double*)malloc((n_row*n_row)* sizeof(double));
  if ((matA == NULL) || (matB == NULL) || (matC == NULL)) {
    exit(1);
  }
  // make arbitrary matrices A and B
  for (i=0;i<n_row;i++) {
    for (j=0;j<n_row;j++) {
      matA[i + n_row*j]=i/(j+3.5);
      matB[i + n_row*j]=j/(i+2.7);
    }
  }
  // initialize timing elements
  double flopings=2.0*n_row*n_row*n_row;
  double tod_start, t_cpu, real_time;
  long start_point, stop_point,total_time;
  float gflops;
  // perform AB r times
  for(l=0;l<n_rep;l++) {
    start_point=readTSC();
    tod_start= gtod();
    //use dgemm_ to calculate the matrix multiplication:
    int N=n_row;
    char TRANSA = 'N';
    char TRANSB = 'N';
    int M=N;
    int K=N;
    double ALPHA=1.;
    double BETA=0.;
    int LDA=N;
    int LDB=N;
    int LDC=N;
    dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,matA,&LDA,matB,&LDB,&BETA,matC,&LDC);
    // determine and display timing results
    real_time= gtod() - tod_start;
    stop_point= readTSC();
    total_time= stop_point- start_point;
    printf("%ld,%f\n",total_time,real_time);
  }
  //free all the memory:
  free(matA);
  free(matB);
  free(matC);

  return 0;
}
