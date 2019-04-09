#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"
#include "cblas.h"
//reference:
// 1. http://www.netlib.org/clapack/cblas/dgemm.c
// 3. https://www.math.utah.edu/software/lapack/lapack-blas/dgemm.html
// 4  http://www.math.utah.edu/software/lapack/lapack-d/dpotrf.html
// 6. https://linux.die.net/man/3/getline
// 7. http://www.cplusplus.com/reference/cstdio/snprintf/
// 8. http://www.netlib.org/lapack/explore-html/d1/d54/group__double__blas__level3_ga6a0a7704f4a747562c1bd9487e89795c.html#ga6a0a7704f4a747562c1bd9487e89795c
// 9. http://gnuplot.sourceforge.net/
// 10.https://stackoverflow.com/questions/30762001/c-gnuplot-pipe-input-from-c-defined-variables

int main(int argc,char *argv[]){
    //Global Constants:
    const char data_directory[] = "./data.dat";
    int d,i,j,n,nrow; 
    char string[1024];
    double x, y,mean,SSE,SSTO,R2;//R2 are the R squared statistic for polynomial fitting.
    FILE *ifp; 
    FILE *raw;
    FILE *fit;
    FILE *gnuplotPipe;
    //The inputs of the LAPACK/BLAS package:
    char TRANSA;
    char TRANSB;
    char TRANS;
    double ALPHA;
    double BETA;
    int M;
    int N;
    int K;
    int LDA;
    int LDB;
    int LDC;
    int INCB;
    int INCY;
    int INCP;
    char SIDE;
    char DIAG;
    /////
	//read in the degree of the polynomial fitted
    d=atoi(argv[1])+1;
    ifp=fopen(data_directory,"r");
    if(ifp==NULL){
        exit(1);
    }
    int degree=d;
    //read the number of lines in the data:
    fscanf(ifp,"%d",&n);
    //xx store the x datalines;
    nrow=n;
    double *xx=(double*)malloc(nrow*sizeof(double));
    //X store the design matrix:
    double *X=(double*)malloc(nrow*d*sizeof(double));
    //Y
    double *Y=(double*)malloc(nrow*sizeof(double));
    for(i=0;i<nrow;i++){
        fscanf(ifp,"%lf,%lf",&xx[i],&Y[i]);
        for(j=0;j<d+1;j++){
            X[i+j*n]=pow(xx[i],j);
        }
    }
    //Calculate A=X^{T}X by dgemm function 
    TRANSA='T'; //X^{T} transform!
    TRANSB='N';
    ALPHA=1.0;
    BETA=0.0;
    M=d;
    N=d;
    K=n; 
    LDA=n;
    LDB=n; 
    LDC=d;
    double *A=(double*) malloc(d*d* sizeof(double));
    if (A == NULL) {
    	exit(1);
    }
    dgemm_(&TRANSA,&TRANSB,&M,&N,&K,&ALPHA,X,&LDA,X,&LDB,&BETA,A,&LDC);



    //Calculate P=X'Y using dgemv:
    TRANS='T';
    M=n;
    N=d;
    LDA=n;
    INCY=1;
    INCP=1;
    double *P=(double*) malloc(d* sizeof(double));
    if(P==NULL){
    	exit(1);
    }
    dgemv_(&TRANS,&M,&N,&ALPHA,X,&LDA,Y,&INCY,&BETA,P,&INCP);


    //Calculate the Choelsky decomposition XTX = A = LL' i.e. get L
    double *L=(double*)malloc(d*d*sizeof(double));
    if(L==NULL){
        exit(1);
    }
    //set the initial value of L
    for(i=0;i<d*d;i++){
        L[i] = A[i];
    }
    char UPLO='L'; // A=LL'
    int ORD=d;//order
    int LDD=d;//dimension
    int INFO=0;
    dpotrf_(&UPLO,&ORD,L,&LDD,&INFO);



    //INFO!=0 means the decomposition failed
    if(INFO!=0){
        exit(1);
    }
    //Solving the linear system Ab=P is indeed a two step problem
    //Sub-step 1:solve LQ=P where L is the lower-triangular cholesky decomposition of A=X^{T}X 
    //Q=L^{-1}P using dtrsm in LAPACK to solve it
    SIDE='L'; //left
    UPLO='L'; // lower triangular same as in Cholesky decomposition
    TRANSA='N'; 
    DIAG='N'; // not diagnoal
    M=d;
    N=1;
    ALPHA=1.0;
    LDA=d;
    LDB=d;
    //dynamic allocate memory for matrix Q:
    double *Q=(double*) malloc(d* sizeof(double));
    if(Q==NULL){
    	exit(1);
    }
    //Set the initial value of Q as P
    for(i=0;i<d;i++){
        Q[i]=P[i];
    }
    dtrsm_(&SIDE,&UPLO,&TRANSA,&DIAG,&M,&N,&ALPHA,L,&LDA,Q,&LDB);


    //Sub-step 2:solve L'B=Q where b is the result for least square problem:
    //B=L'^{-1}Q
    SIDE='L';
    UPLO='L'; 
    TRANSA='T';
    DIAG='N';
    M=d;
    N=1;
    ALPHA=1.0; 
    LDA=d;
    LDB=d;
    double *B=(double*) malloc(d* sizeof(double));
    if (B==NULL) {
    	exit(1);
    }
    //Set the initial value of B
    for(i=0;i<d; i++){ 
        B[i]=Q[i];
    }
    dtrsm_(&SIDE,&UPLO,&TRANSA,&DIAG,&M,&N,&ALPHA,L,&LDA,B,&LDB);

    //Calculate the fitting value Yhat=XB:
    double *Yhat = (double*) malloc(n* sizeof(double));
    if (Yhat == NULL) {
    	exit(1);
    }
    TRANS='N';
    M=n;
    N=d;
    LDA=n;
    ALPHA=1.0;
    BETA=0.0;
    INCB=1;
    INCY=1;
    dgemv_(&TRANS,&M,&N,&ALPHA,X,&LDA,B,&INCB,&BETA,Yhat,&INCY);

    /////////////////////////
    //print all the data we get through matrix computing::
    //print the matrix A=X'X:
    printf("PRINT MATRIX A=XtX:\n");
    for(i=0;i<d;i++){
        for(j=0;j<d;j++){
            printf("%lf\t",A[i+j*d]);
        }
        printf("\n");
    }
    printf("\n");
    //print the matrix P:
    printf("PRINT MATRIX P:\n");
    for(i=0;i<d;i++){
        printf("%lf\n",P[i]);
    }
    printf("\n"); 
    //print the matrix L:
    printf("PRINT MATRIX L:\n");
    for(i=0;i<d;i++){
        for(j=0;j<d;j++){
            if(i<j){
                L[i+j*d]=0.00;
            }
            printf("%lf\t",L[i+j*d]);
        }
        printf("\n");
    }
    printf("\n");
    //print the matrix Q:
    printf("PRINT MATRIX Q:\n");
    for(i=0;i<d;i++){
        printf("%lf\n",Q[i]);
    }
    printf("\n");
    //print the matrix B:
    printf("PRINT MATRIX B:\n");
    for(i=0;i<d;i++){
        printf("%lf\n",B[i]);
    }
    printf("\n");    
    //print the matrix Yhat:
    printf("PRINT MATRIX YHAT:\n");
    for(i=0;i<n;i++){
        printf("%lf\n",Yhat[i]);
    }



    //Calculate R2 statistic to measure the fitting:
    mean=0;
    SSE=0.0;
    SSTO=0.0;
    for(i=0;i<n;i++){
    	mean+=Y[i]/n;
    	SSE+=pow((Yhat[i]-Y[i]),2);
    }
    for(i=0;i<n;i++){
    	SSTO+=pow((Y[i]-mean),2);
    }
    R2=1.0-SSE/SSTO;
    //output the results as files for downstream plotting:
    snprintf(string,sizeof(string),"%d%s",degree-1,"poly_design.dat");
    raw=fopen(string,"w");// raw store the Y data points values~ X point values
    snprintf(string,sizeof(string),"%d%s",degree-1,"poly_fitting.dat");
    fit=fopen(string,"w");// fit store the fitted Y values~ X point values
    for (i=0; i < n; i++) {
        X[i+n]=X[i+n];
        Y[i]=Y[i];
        fprintf(raw,"%lf %lf\n",X[i+n],Y[i]);
        fprintf(fit,"%lf %lf\n",X[i+n],Yhat[i]);
    }
    //close the files
    fclose(raw);
    fclose(fit);

    //// Gnuplot the resulting fit and the raw data in one plot:
    gnuplotPipe = popen("gnuplot", "w");//call gnuplot interface in C
    fprintf(gnuplotPipe,"set terminal jpeg\n");//output the jpeg picture
    fprintf(gnuplotPipe,"set output './plot_degree_Heqiao%d.png'\n", d-1); //set the name of plot
    fprintf(gnuplotPipe,"set title 'Polynomial Fitting points and observed point with degree %d\n",d-1);// set title 
    fprintf(gnuplotPipe,"set xlabel 'X Datapoints'\n" );
    fprintf(gnuplotPipe,"set ylabel 'Y Datapoints'\n" );//set X and Y label
    fprintf(gnuplotPipe,"set key left top\n");//set the legend at the top left
    fprintf(gnuplotPipe,"set pointsize 3\n" );//set the point size
    fprintf(gnuplotPipe,"plot './%dpoly_design.dat' title 'Input', ", d-1);//plot the raw (X_i,Y_i) data points
    fprintf(gnuplotPipe,"'./%dpoly_fitting.dat' lt rgb 'violet' title 'Fit'\n", d-1);//plot the fitted (X_i,Y_i) data points
    pclose(gnuplotPipe);    
    return 0;
}
