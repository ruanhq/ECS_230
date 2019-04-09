//Compute the leading eigenvector:
#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "cblas.h"
//Reference:
// 1. http://www.netlib.org/clapack/cblas/dgemm.c
// 2. https://www.math.utah.edu/software/lapack/lapack-blas/dgemm.html
// 3.  http://www.math.utah.edu/software/lapack/lapack-d/dpotrf.html
// 4. https://linux.die.net/man/3/getline
// 5. http://www.cplusplus.com/reference/cstdio/snprintf/
// 6. http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_ga6a0a7704f4a747562c1bd9487e89795c.html#ga6a0a7704f4a747562c1bd9487e89795c
// 7. http://gnuplot.sourceforge.net/
// 8.https://stackoverflow.com/questions/30762001/c-gnuplot-pipe-input-from-c-defined-variables
// 9.http://www.netlib.org/lapack/explore-html/de/da4/group__double__blas__level1_ga4b25c539b862414e6f91ebb30b391d7c.html#ga4b25c539b862414e6f91ebb30b391d7c
// 10.https://stackoverflow.com/questions/3118165/when-do-i-use-fabs-and-when-is-it-sufficient-to-use-stdabs


const char data_directory[] = "./figuring.txt";
int i,j,k,itt;
int n;
int line_1;
FILE *ifp;
char TRANS;
double ALPHA;
double BETA;
int M;
int N;
int LDA;
int LDB;
int LDC;
int INCA;
int INCB;
int main(int argc, char** argv)
{
    // read data into memory from disk
    ifp = fopen(data_directory, "r");
    if (ifp == NULL) {
        printf("OPEN FILE FAILED\n"); // error if it can't open
        exit(1);
    }
    n = 1;
    line_1=1;
    //memory allocation:
    double *A = (double*) malloc(n*n* sizeof(double)); 
    double *b1 = (double*) malloc(n* sizeof(double)); 
    double *b2 = (double*) malloc(n* sizeof(double)); 
    // read file
    char * line = NULL;
    size_t len = 0;
    ssize_t read;
    char * pch;
    i = 0;
    while ((read = getline(&line, &len, ifp)) != -1){
        if(line_1== 1){
            n = atoi(line); // n_obs is first line of data.dat
            line_1= 0;
            /* printf("Using %d observations\n", n); */
            printf("Observation n=%d\n", n);
            // allocate memory for the input data
            // NOTE: A is col-major indexed i.e. A[i + n*j] = A_(i,j)
            A = realloc(A, n*n*sizeof(double));
            b1 = realloc(b1, n*sizeof(double));
            b2 = realloc(b2, n*sizeof(double));
            if (A == NULL) {
                printf("ALLOCATION FAIL!\n");
                return(1);
            }else{
            }
            int i1, i2;
            for(i1=0;i1<n;i1++){
                for(i2=0;i2<n;i2++){
                    A[i1 + n*i2]=0.0;
                }
                b1[i1] = (double)i1+1;
                b2[i1] = b1[i1];
            }
        }else{
            pch = strtok(line," ");//split the string by the seperation " "
            k = 0;
            line_1= 1;
            double nj;
            while (pch != NULL) {
                if(line_1== 1){
                    nj = atof(pch); 
                    line_1= 0;
                }else{
                    k = atof(pch);
                    A[k-1 + n*i] = 1.0/nj;//calculate it as the reciprocal of the degree
                }
                pch = strtok(NULL, " ");
            }
            i++;
            /* printf("\n"); */
        }
    }
    fclose(ifp);
    printf("\n");
    i=0;
    j=0;
    // print A
    printf("PRINT MATRIX:A:\n");
    for(i=0;i<n;i++){
        for(j=0;j<n;j++){
                printf("%f ", A[i + n*j]);
        }
        printf("\n");
    }
    //shift:
    //double translation=0.3;
    //for(i=0;i<n;i++){
    //    A[i+n*i]=A[i+n*i]+translation;
    //}
    //modified version
    //double *S=(double*)malloc(n*n*sizeof(double));
    //double *MM=(double*)malloc(n*n*sizeof(double));
    //for(i=0;i<n*n;i++){
    //    for(j=0;j<n;j++){
    //        S[i+n*j]=(double)1.0/n;
    //    }
    //}
    //for(i=0;i<n;i++){
    //    for(j=0;j<n;j++){
    //        MM[i+n*j]=0.85*A[i+n*j]+0.15*S[i+n*j];
    //    }
    //}
    //print matrix MM:
    //printf("PRINT MATRIX MM:\n");
    //for(i=0;i<n;i++){
    //    for(j=0;j<n;j++){
    //        printf("%f",MM[i+n*j]);
    //    }
    //    printf("\n");
    //}
    //normalize the b:
    //initialize the b as [1,1,1,1]
    double norm_b;
    int INCX=1;
    norm_b=cblas_dnrm2(n,b1,INCX);
    for(i=0;i<n;i++){
        b1[i]=b1[i]/norm_b;
    }
    norm_b=cblas_dnrm2(n,b1,INCX);
    //print the normalized one:
    printf("\n PRINT Normalized b0:\n");
    for(i=0;i<n;i++){
        printf("%f\n",b1[i]);
    }
    printf("\n");
    //
    //compute:
    //specify the constant to use:
    TRANS = 'N';
    M = n; 
    N = n;
    LDA = n;
    INCB = 1;
    ALPHA = 1.0;
    BETA = 0.0;
    double epsilon = 1.0e-6;
    double delta = 1.0;
    int breaker = 150;
    itt = 0;
    while(delta > epsilon){
        // copy b1 <- b2
        for(i=0;i<n;i++){
            b1[i] = b2[i];
        }
        // normalize b1
        norm_b = cblas_dnrm2(n, b1, INCX);
        for(i=0; i<n; i++){
            b1[i] = b1[i] / norm_b;
        }
        // multiply b2 <- A b1
        dgemv_(&TRANS, &M, &N, &ALPHA,MM, &LDA,b1, &INCB, &BETA,b2, &INCB);
        // normalize b2
        norm_b = cblas_dnrm2(n, b2, INCX);
        for(i=0;i<n;i++){
            b2[i] = b2[i] / norm_b;
        }
        // calculate difference
        delta = 0.0;
        for(i=0; i<n; i++){
            delta = delta + fabs(b1[i] - b2[i]);
        }
        printf("Iteration Number %d: the difference is %f\n", itt, delta);

        itt++;
        if(itt > breaker){
            printf("\n Still not converge!\n");
            break;
        }
    }
    //Here we scale the vector so that the norm becomes 1.
    norm_b = cblas_dasum(n, b1, 1);
    for(i=0; i<n; i++){
        b1[i] = b1[i] / norm_b;
    }


    if(itt>breaker){
        printf("Can't converge After %d iterations to vectors:\n", itt);
        printf("The vector now is:\n");
        for(i=0; i<n; i++){
            printf("%f\n", b2[i]);
        }
    }else{
        printf("Converge after %d iterations,the eigenvector we solve here is:\n", itt);
        for(i=0; i<n; i++){
            printf("%f\n", b1[i]);
        }
    }
    return 0;
}



