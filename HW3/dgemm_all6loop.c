#include "stdio.h"
#include "time.h"
#include "sys/time.h"
#include "math.h"
#include "cblas.h"
#include "stdlib.h"

//Reference:
//1. https://canvas.ucdavis.edu/courses/246067/files/folder/Homework/timing1.c
//2. https://canvas.ucdavis.edu/courses/246067/files/folder/Homework/timing2.c
//3. https://stackoverflow.com/questions/39002052/how-i-can-print-to-stderr-in-c
//4. https://stackoverflow.com/questions/328834/c-delete-vs-free-and-performance
//5. https://stackoverflow.com/questions/39002052/how-i-can-print-to-stderr-in-c
//6. https://en.wikipedia.org/wiki/C_dynamic_memory_allocation
//7. https://stackoverflow.com/questions/28654438/matrix-vector-product-with-dgemm-dgemv
//8. https://www.tutorialspoint.com/c_standard_library/c_function_clock.htm
//9. http://www.netlib.org/blas

//cite from canvas code:
long long readTSC(void)
{
 /* from timing2.c downloaded in canvas */
  union { long long complete; unsigned int part[2]; } ticks;
  __asm__ ("rdtsc; mov %%eax,%0;mov %%edx,%1"
            : "=mr" (ticks.part[0]),
              "=mr" (ticks.part[1])
            : /* no inputs */
            : "eax", "edx");
  return ticks.complete;
}

//get the time of the day directly downloaded from Professor Gygi's sample code
volatile double gtod(void) {
  static struct timeval tv;
  static struct timezone tz;
  gettimeofday(&tv,&tz);
  return tv.tv_sec + 1.e-6*tv.tv_usec;
}



//
int main(int argc, char** argv){
	if((argc<2)||(argc>3)){
		return(-1);
	}
	clock_t clock_start,clock_end;
	double tod_start,tod_dif;
	long starting_point,stoping_point,total_time;
	float gflops;
	int n=atoi(argv[1]);
	int r=atoi(argv[2]);
	int i,j,k,l;
	double sum=0.0;
	//dynamic memory allocation:
	double *matA=(double*)malloc((n*n)*sizeof(double));
	double *matB=(double*)malloc((n*n)*sizeof(double));
	double *matC=(double*)malloc((n*n)*sizeof(double));
	if((matA==NULL)||(matB==NULL)||(matC==NULL)){
		return(-1);
	}
	long flopings=2*n*n*n;
	//random initialize:
	for(i=0;i<n;i++){
		for(j=0;j<n;j++){
			matA[i+n*j]=j/(i+3.05)+3.48/(i+n*j+0.98);
			matB[i+n*j]=i/(j+1.05)+0.05/(i+n*j+0.002);
		}
	}	
	for(l=0;l<r;l++){
		clock_start=clock();
		starting_point=readTSC();
		tod_start=gtod();
		//possible order 1:jik
		for(j=0;j<n;j++){
			for(i=0;i<n;i++){
				for(k=0;k<n;k++){
					sum=sum+matA[i+n*k]*matB[k+n*j];
				}
				matC[i+n*j]=sum;
				sum=0.0;
			}
		}
		clock_end=clock();
		long long clock_dif1=clock_end-clock_start;
		tod_dif=gtod()-tod_start;
		stoping_point=readTSC();
		total_time=stoping_point-starting_point;
		gflops=0.000000001*flopings/tod_dif;
		printf("%ld,%f,%f possible order1:jik \n",total_time,tod_dif,gflops);
		//possible order 2:ikj:
		clock_start=clock();
		starting_point=readTSC();
		tod_start=gtod();
		for(j=0;j<n;j++){
			for(k=0;k<n;k++){
				for(i=0;i<n;i++){
					sum=sum+matA[i+n*k]*matB[k+n*j];
				}
				matC[i+n*j]=sum;
				sum=0.0;
			}
		}
		clock_end=clock();
		long long clock_dif2=clock_end-clock_start;
		tod_dif=gtod()-tod_start;
		stoping_point=readTSC();
		total_time=stoping_point-starting_point;
		gflops=0.000000001*flopings/tod_dif;
		printf("%ld,%f,%f possible order2:jki \n",total_time,tod_dif,gflops);
		//jki:
		clock_start=clock();
		starting_point=readTSC();
		tod_start=gtod();
		for(i=0;i<n;i++){
			for(j=0;j<n;j++){
				for(k=0;k<n;k++){
					sum=sum+matA[i+n*k]*matB[k+n*j];
				}
				matC[i+n*j]=sum;
				sum=0.0;
			}
		}
		clock_end=clock();
		long long clock_dif3=clock_end-clock_start;
		tod_dif=gtod()-tod_start;
		stoping_point=readTSC();
		total_time=stoping_point-starting_point;
		gflops=0.000000001*flopings/tod_dif;
		printf("%ld,%f,%f, possible order3:ijk \n",total_time,tod_dif,gflops);
		//jik:
		clock_start=clock();
		starting_point=readTSC();
		tod_start=gtod();
		for(k=0;k<n;k++){
			for(j=0;j<n;j++){
				for(i=0;i<n;i++){
					sum=sum+matA[i+n*k]*matB[k+n*j];
				}
				matC[i+n*j]=sum;
				sum=0.0;
			}
		}
		clock_end=clock();
		long long clock_dif4=clock_end-clock_start;
		tod_dif=gtod()-tod_start;
		stoping_point=readTSC();
		total_time=stoping_point-starting_point;
		gflops=0.000000001*flopings/tod_dif;
		printf("%ld,%f,%f, possible order4:kji \n",total_time,tod_dif,gflops);
		//kij:
		clock_start=clock();
		starting_point=readTSC();
		tod_start=gtod();
		for(i=0;i<n;i++){
			for(k=0;k<n;k++){
				for(j=0;j<n;j++){
					sum=sum+matA[i+n*k]*matB[k+n*j];
				}
				matC[i+n*j]=sum;
				sum=0.0;
			}
		}
		clock_end=clock();
		long long clock_dif5=clock_end-clock_start;
		tod_dif=gtod()-tod_start;
		stoping_point=readTSC();
		total_time=stoping_point-starting_point;
		gflops=0.000000001*flopings/tod_dif;
		printf("%ld,%f,%f possible order5:ikj \n",total_time,tod_dif,gflops);
	    //kji:
		clock_start=clock();
		starting_point=readTSC();
		tod_start=gtod();
		for(k=0;k<n;k++){
			for(i=0;i<n;i++){
				for(j=0;j<n;j++){
					sum=sum+matA[i+n*k]*matB[k+n*j];
				}
				matC[i+n*j]=sum;
				sum=0.0;
			}
		}
		clock_end=clock();
		long long clock_dif6=clock_end-clock_start;
		tod_dif=gtod()-tod_start;
		stoping_point=readTSC();
		total_time=stoping_point-starting_point;
		gflops=0.000000001*flopings/tod_dif;
		printf("%ld,%f,%f possible order6:kij \n",total_time,tod_dif,gflops);    
	}
	//remove all the memory:
	free(matA);
	free(matB);
	free(matC);
	return 0;
}














