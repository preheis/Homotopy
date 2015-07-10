#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//#include <cblas.h>
//#include <lapacke.h>
#include <time.h>

/*This program computes all of the eigenpairs of symmetric tridiagonal matricies using a Homotopy algorithm*/

//Global variables
int n = 100;
int k = 0;

/*This function reduces the step size, alters the homotopy parameter, and resets the eigenvector*/
int ReduceStep(double *A, double *D, double *DD, double *XT, double *DDD, double *OO, double H, double PT, double Z, double Z1, double APP, int NS){
	double *X = (double*)calloc(n,sizeof(double));
	double T;
	H = H/2;
	T = PT+H;
	cblas_dcopy(n,XT,1,X,1);
	if (NS==0){
		double PE = Z+Z1*H;
		mainblock(A, D, DD, DDD, OO, X, XT, Z, Z1, H, T, PT, PE, NS);
	}
	else{
		Predict(A, D, DD, X, DDD, OO, H, T, Z, Z1, APP, NS);
	}		
	return 0;	
}

/*Function recomputes PE using Hermite Interpolation.*/
/*If T == 1, the eigenvalues are stored.*/
int Predict(double *A, double *D, double *DD, double *X, double *DDD, double *OO, double H, double T, double Z, double Z1, double APP, int NS){
	double *XT = (double*)calloc((double)n,sizeof(double));
	double *HA = (double*)calloc((double)n,sizeof(double));
	double *EV = (double*)calloc((double)n,sizeof(double));
	double *EVT= (double*)calloc((double)n*n,sizeof(double));	
	double OZ=Z, OZ1=Z1, Q, QQ, PH=H, PT=T,PE;
	int i;
	H = 1-PT; T = 1; NS+=1; 
	cblas_dcopy(n,X,1,XT,1);
	if (T==1){
		EV[k]=Z;
			
	}
	else{
		//Recompute Z1
		HA[0]=DDD[0]*X[0]+OO[1]*X[1];
		for (i = 1; i<n-1; i++){
			HA[i]=OO[i]*X[i-1]+DDD[i]*X[i]+OO[i+1]*X[i+1];
		}
		HA[n-1]=OO[n-1]*X[n-2]+DDD[n-1]*X[n-1];
		//Recompute PE
		Q = pow(1+(PH/H),2);
		QQ = Q*(H/PH);
		PE=OZ+OZ1*(H+PH)+Q*((Z-OZ)-OZ1*PH)+QQ*(PH*(Z+OZ1)-2*(Z-OZ));

		mainblock(A, D, DD, DDD, OO, X, XT, Z, Z1, H, T, PT, PE, NS);
	}
	return 0;
}

/*Function generates the type of matrix as described in Example 1*/
/*int Example1(double *A, int n, int m, int k){
	int i,j,r;
	double l;
	for (i = 0; i<n; i++){
		A[i*n+i]= (double) i;
		for (j = 1; j<n; j++){
			r = j%2;
			if (r==0){
				l = sqrt(j-(0.5*j));
				A[j*n+(j-1)]=l;
				//A[(j-1)+n*j]=l;
			}
			else{
			    l = k*sqrt(m+(j-1)-(0.5*(j-1)));
				A[j*n+(j-1)]=l;
				//A[(j-1)+n*j]=l;
			}						
		}
	}
	return 0;
}*/

/*Function generates the type of matrix as described in Example 3*/
int Example3(double *A){
	int i,j;
	for(i = 0; i<n; i++){
		A[i*n+i]= (double) i;
		for (j = 1; j<n; j++){
			A[j*n+(j-1)]=1.0;
			//A[(j-1)+j*n]=1.0;		
		}	
	}
	return 0;
}

/*Function generates an identity matrix of size n*/
int IdentityMat(double *ID){
	int i;
	double j=1;
	for (i = 0; i<n; i++){
		ID[(i*n)+i] = j;	
	}
	return 0;
}


/*Initialize all of the vectors needed for computation*/
int Initialize(double *A, double *D, double *DD, double *DDD, double *OO){
	int i;
	for (i = 0; i<n; i++){
		DD[i] = A[i*n+i];
		D[i*n+i] = DD[i];
		DDD[i] = DD[i]-DD[i];
	}
	OO[0]=0;
	for (i = 1; i<n; i++){
		OO[i]=A[i*n+(i-1)];
	}
	return 0;
}

/*Compute the number of sign changes using Sturm Sequence method*/
int COUNT(double *DD, double *OO, double x){
	int Count = 0, d = 1,i;
	
	for (i = 0; i<n; i++){
		double b = pow(OO[i],2);		
		d = (DD[i]-x-b)/d;
		if (d<0){
			Count+=1;
		}		
	}
	return Count;
}

/*Approximate an eigenvector with Inverse Iteration*/
double * II(double *At, double *X, double APP, int j){
	double *W = (double*)calloc(n*n,sizeof(double));
	double *ID = (double*)calloc(n*n,sizeof(double));
	double *Y = (double*)calloc(n,sizeof(double));
	int *IPIV = (int*)calloc(n,sizeof(int));
	double norm;	
	int i,NRHS=1,LDY=1, INFO;
	
	IdentityMat(ID);
    cblas_dcopy(n,At,1,W,1); //W = At;
	cblas_dscal(n,APP,ID,1); //I = I*APP;
	cblas_daxpy(n,-1,ID,1,W,1); // W = W-I;
	
	for (i=0; i<=j; i++){
		//Perform Inverse Iteration
		norm = cblas_dnrm2(n,X,1); //norm(X);
		cblas_dscal(n,1/norm,X,1); // X <- X*1/norm(X)
		cblas_dcopy(n,X,1,Y,1); // Y <- X
	    INFO = dgesv_(n,NRHS,W,n,IPIV,Y,LDY); //W*X = Y 		
	}
	norm = cblas_dnrm2(n,X,1); //norm(X);
	cblas_dscal(n,1/norm,X,1); // X <- X*1/norm(X)
    cblas_dcopy(n,X,1,Y,1); // Y <- X
	free(W);
	free(ID);
	return Y;
}

/*Approximate an eigenvalue with Rayleigh Quotient*/
double RQI(double *At, double *X, int j){
	double *u = (double*)calloc(n,sizeof(double));
	double norm,APP;
	int i;
	int incx=1;
	int alpha=1;
	int beta = 1;
	int LDA=1;
	int incy=1; 
	char N;
	for (i = 0; i<=j; i++){
		//Peform Rayleigh Quotient Iteration
		norm = cblas_dnrm2(n,X,incx); //norm(X)
		cblas_dscal(n,1/norm,X,incx); //X <- X*1/norm(X);
		dgemv_(N,n,n,alpha,*At,n,*X,incx,beta,u,&incy); //u <- At*X
		APP = cblas_ddot(n,X,incx,u,incy); //APP <- X'*u
	}
	free(u);
	return APP;
}
/*The main block of the program. It corrects the prediction computed.*/
int mainblock(double *A, double *D, double *DD, double *DDD, double *OO, double *X, double *XT, double Z, double Z1, double H, double T, double PT, double PE, int NS){

	double eps, EPS, EPS1, EPS2, EPS3, NORM, SUMA, SUMB, APP, RES;
	int i,sc,KK;
	
	//Set the error tolerances for the program
	eps = 1.11E-16;

	for (i = 0; i<n; i++){
		SUMA+=abs(DD[i]);
		SUMB+=abs(OO[i]);
	}
	
	NORM = SUMA+SUMB;
	
	EPS=eps*n*NORM;
	EPS2=EPS3=EPS;
	if (T==1){
		EPS1=EPS;
	}
	else{
		EPS1=sqrt(EPS);
	}
	
	//////////////////////////////////////////////////

	double *At = (double*)calloc(n*n,sizeof(double));
	double *Y = (double*)calloc(n,sizeof(double));
	double *AY = (double*)calloc(n,sizeof(double));
	double *APPY = (double*)calloc(n,sizeof(double));	
	int j=1;
	int incx=1;
	int incy=1;
	int nincx=-1;
	int alpha =1;
	int beta = 1;
	int LDA=1;
	char N;
	//At = D + T*(A-D);
	cblas_dcopy(n,A,1,At,1);//At <- A
	cblas_daxpy(n,nincx,D,incx,At,incy); //At <- A-D
	cblas_dscal(n,T,At,incx); //At <- T*(A-D)
	cblas_daxpy(n,alpha,D,incx,At,incy); // At <- D + T*(A-D)
	cblas_dcopy(n,X,1,Y,1);
	APP = PE;
	sc = COUNT(DD,OO,PE);
	
	if (sc!=k && sc!=k-1){ReduceStep(A, D, DD, XT, DDD, OO, H, PT, Z, Z1, APP, NS); return 0; }	
	if (sc==k){KK=1;}
	if (sc==k-1){KK=0;}
	
	for (i = 0; i<=10; i++){
		Y = II(At,Y,APP,j);
		//Compute Residual
		cblas_dcopy(N,Y,incy,APPY,incy);//APPY <- Y
		dgemv_(N,n,n,alpha,A,LDA,Y,incy,beta,AY,incy);//AY <- A*Y
		cblas_dscal(n,APP,APPY,incy);//APPY <- APP*Y
		cblas_daxpy(N,-alpha,APPY,incy,AY,incy);//RES <- A*Y-APP*Y
		RES = isamax_(n,APPY,1);//RES <- max(RES);
		if (RES<=EPS1){
			sc = COUNT(DD,OO,APP-EPS1-EPS3);
			if (sc==k-1){
				Predict(A, D, DD, X, DDD, OO, H, T, Z, Z1, APP, NS);
				return 0;			
			}
			else if (sc>k-1){
				ReduceStep(A, D, DD, XT, DDD, OO, H, PT, Z, Z1, APP, NS);
				return 0;
			}
			else{
				sc = COUNT(DD,OO,APP+EPS1+EPS3);
				if (sc>=k){
					Predict(A, D, DD, X, DDD, OO, H, T, Z, Z1, APP, NS);
					return 0;				
				}
			}		
		}
		
		APP = RQI(At,Y,j);
		
		if (KK==1 && APP>(PE+EPS2)){
			ReduceStep(A, D, DD, XT, DDD, OO, H, PT, Z, Z1, APP, NS);
			return 0;
		}
		
		if (KK==0 && APP<(PE-EPS2)){
			ReduceStep(A, D, DD, XT, DDD, OO, H, PT, Z, Z1, APP, NS);
			return 0;
		}
		
		if (i==10){
			ReduceStep(A, D, DD, XT, DDD, OO, H, PT, Z, Z1, APP, NS);
			return 0;
		}
	}
	free(At);
	free(Y);
	free(AY);
	free(APPY);
	return 0;		
}
//Main
int main(void){
	//Other important variables
	double alpha, beta, H, PT, T, Z, Z1, Z2, Z3, PE;	
	int NS=0,i,incx, incy,nincx;

	H = PT = 0;
	T = 1.0;

	incx = incy = 1;
	nincx=-1;
	alpha = beta = 1.0;
	clock_t t; 
	//Matrices
	double *A = (double*)calloc(n*n,sizeof(double));
	double *D = (double*)calloc(n*n,sizeof(double));
	double *Acopy = (double*)calloc(n*n,sizeof(double));
	//Vectors	
	double *d = (double*)calloc(n,sizeof(double));
	double *DD = (double*)calloc(n,sizeof(double));
	double *DDD = (double*)calloc(n,sizeof(double));
	double *OO = (double*)calloc(n,sizeof(double));
	double *X = (double*)calloc(n,sizeof(double));	
	double *XT = (double*)calloc(n,sizeof(double));
	double *XTprime = (double*)calloc(n,sizeof(double));

	Example3(A);
	Initialize(A,D,DD,DDD,OO);
	d[k]=1.0;
	t = clock();
	printf("Beginning Computation\n");
	for (i = 0; i<n; i++){
		k=i;
		Z=d[k];
		XT[k]=1;
		printf("Copying Vectors and Matrices\n");
		//Copy vectors and matrices
	    cblas_dcopy(n,A,incx,Acopy,incy);// Acopy <- A
		cblas_dcopy(n,XT,incx,X,incx);//X <- XT
		//Compute first derivative
		printf("Computing Z1\n");
		cblas_daxpy(n,nincx,D,incx,Acopy,incx);//A <- -1*D + A
		printf("!!!!\n");
		dgemv_("t",&n,&n,&alpha,Acopy,&n,XT,&incx,&beta,XTprime,&incy); //XT' <- (A-D)*XT
		printf("!!!!\n");
		Z1 = cblas_ddot(n,XT,1,XTprime,1);//Z1 <- XT'*((A-D)*XT)
		//Predict eigenvalue with first order taylor method
		PE=Z+Z1*H;
		printf("!!!!\n");
		printf("PE = %3.5f\n",PE);
		printf("Beginning Correction.\n");
		mainblock(A, D, DD, DDD, OO, X, XT, Z, Z1, H, T, PT, PE, NS);
	}
	//Free all allocated memory
	free(A);
	free(D);
	free(d);
	free(DD);
	free(DDD);
	free(OO);
	free(X);
	free(XT);
	t = clock()-t;
	//fprintf("Total execution time: %f seconds\n", (float)t/CLOCKS_PER_SEC);
	return 0;
}


