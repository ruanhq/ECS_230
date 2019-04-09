#include "stdio.h"
#include "time.h"
#include "sys/time.h"
#include "math.h"
#include "stdlib.h"

//Reference:
//1. https://stackoverflow.com/questions/39002052/how-i-can-print-to-stderr-in-c
//2. https://en.wikipedia.org/wiki/C_dynamic_memory_allocation
//3. https://stackoverflow.com/questions/28654438/matrix-vector-product-with-dgemm-dgemv
//4. https://www.tutorialspoint.com/c_standard_library/c_function_clock.htm
//5. https://canvas.ucdavis.edu/courses/246067/files/folder/Homework/timing1.c
//6. https://canvas.ucdavis.edu/courses/246067/files/folder/Homework/timing2.c
//7. http://www.cplusplus.com/reference/cstdlib/atoi/

//readTSC function to read the time stamp counter on the chips:

long long readTSC(void)
{
 /* read the time stamp counter on Intel x86 chips */
 /* from timing2.c downloaded in canvas */
  union { long long complete; unsigned int part[2]; } ticks;
  __asm__ ("rdtsc; mov %%eax,%0;mov %%edx,%1"
            : "=mr" (ticks.part[0]),
              "=mr" (ticks.part[1])
            : /* no inputs */
            : "eax", "edx");
  return ticks.complete;
}


//gtod function from Professor Gygi's sample code.
volatile double gtod(void)
{ 
  /* from timing1.c downloaded in canvas */
  static struct timeval tv;
  static struct timezone tz;
  gettimeofday(&tv,&tz);
  return tv.tv_sec + 1.e-6*tv.tv_usec;
}

//main function to multiply the matrix:
int main(int argc,char** argv){
  //global constants:
  if(argc<2){
    exit(1);
  }
  int i,j,k,l;
  int n_row=atoi(argv[1]);
  int n_col=n_row;
  int n_rep=atoi(argv[2]);
  double sum=0.00;
  //report the time:
  float gflops;
  double tod_start,cpu_time,real_time;
  clock_t clock_start;
  double starting_point;
  long stoping_point,total_time;
  // do the dynamic memory allocation for original ones:
  double *matA=(double*)malloc((n_row*n_col)*sizeof(double));
  double *matB=(double*)malloc((n_row*n_col)*sizeof(double));
  double *matAB=(double*)malloc((n_row*n_col)*sizeof(double));
  if((matA==NULL)||(matB==NULL)||(matAB==NULL)){
    return(-1);
  }
  //randomly initialize the matrices:
  //considering the overflow:
  for(i=0;i<n_row;i++){
    for(j=0;j<n_col;j++){
      matA[i+n_row*j]=j/(i+n_row*j+2.0);
      matB[i+n_row*j]=i/(j+n_row*i+1.0);
    }
  }
  for(l=1;l<n_rep;l++){
    clock_start=clock();
    starting_point=readTSC();
    tod_start=gtod();
    for(i=0;i<n_row;i++){
      for(j=0;j<n_row;j++){
        for(k=0;k<n_row;k++){
          sum=sum+matA[i+n_row*k]*matB[k+n_row*j];
      }
      matAB[i+n_row*j]=sum;
      sum=0.0;
    }
  }
  double flopings=2.0*n_row*n_row*n_row;
  //get the timing:
  long long clock_dif=clock()-clock_start;
  cpu_time=(clock_dif)/CLOCKS_PER_SEC;
  real_time=gtod()-tod_start;
  stoping_point=readTSC();
  total_time=stoping_point-starting_point;
  gflops=0.000000001*flopings/real_time;
  printf("%ld,%f,%f\n",total_time,real_time,gflops);
  }
  return 0;
}




