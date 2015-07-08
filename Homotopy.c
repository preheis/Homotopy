#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cblas.h>
#include <time.h>
/*This program computes all of the eigenpairs of symmetric tridiagonal matricies using a Homotopy algorithm*/

//Global variables
int n=100, k=0;

//Function Declarations
int Example1(double* A, int n);
int Example3(double* A, int n);
int Initialize(double *A, double *D, double* DD, double* DDD, double* d, double* OO, int n);
int mainblock(double* A, double* D, double* DD, double* DDD, double* OO, double* X, double* XT, double Z, double Z1, double H, double T, double PT, double PE, int NS);
int Identity(double* I, int n);
int COUNT(double *At, double* DD, double* OO, double x)
int ReduceStep(double* A, double* D, double* XT, double* DDD, double* OO double H, double PT, double Z, double Z1, double APP);
int Predict(double* A, double* D, double* X, double* DDD, double* OO, double H, double T, double Z, double Z1, double APP);
///////////////////////

//Main
int main(void){
	//Other important variables
	int NS=0,alpha=0;
	clock_t t; 
	double H, PT, T, Z, Z1, Z2, Z3, PE;
	H = PT = 0;
	T = 1;
	//Matrices
	double *A = (double*)calloc(n*n,sizeof(double));
	double *D = (double*)calloc(n*n,sizeof(double));
	double *Acopy = (double*)calloc(n*n,sizeof(double));
	double *Dcopy = (double*)calloc(n*n,sizeof(double));
	//Vectors	
	double *d = (double*)calloc(n,sizeof(double));
	double *DD = (double*)calloc(n,sizeof(double));
	double *DDD = (double*)calloc(n,sizeof(double));
	double *OO = (double*)calloc(n,sizeof(double));
	double *X = (double*)calloc(n,sizeof(double));	
	double *XT = (double*)calloc(n,sizeof(double));
	double *Xcopy = (double*)calloc(n,sizeof(double));
	double *XTcopy = (double*)calloc(n,sizeof(double));

	Example3(A,n);
	Initialize(A,D,DD,DDD,d,OO,n);
	
	t = clock();
	for (int i==0; i<=n; i++){
		k=i;
		Z=d(k);
		XT[k]=1;
		//Copy vectors and matrices
		cblas_dcopy(n,A,1,Acopy,1);
		cblas_dcopy(n,D,1,Dcopy,1);
		cblas_dcopy(n,XT,1,XTcopy,1);
		cblas_dcopy(n,XT,1,X,1);
		cblas_dcopy(n,X,1,Xcopy,1);
		//Perform computations
		cblas_daxpy(n,-1,D,1,A,1);//A <= -1*D + A
		//Multiply result by XT transpose
		//Multiply result by XT


		PE=Z+Z1;
		mainblock(A, D, DD, DDD, OO, X, XT, Z, Z1, H, T, PT, PE, NS);

	}
	free(A);
	free(D);
	free(d);
	free(DD);
	free(DDD);
	free(OO);
	free(X);
	free(XT);
	t = clock()-t;
	fprintf("Total execution time: %f seconds\n", (float)t/CLOCKS_PER_SEC);
	return 0;
}

/*The main block of the program. It corrects the prediction computed.*/
int mainblock(double* A, double* D, double* DD, double* DDD, double* OO, double* X, double* XT, double Z, double Z1, double H, double T, double PT, double PE, int NS){

	double eps, EPS, EPS1, EPS2, EPS3, NORM, SUMA, SUMB, APP, RES;
	int i,sc,KK;
	
	//Set the error tolerances for the program
	eps = 1.11E-16;

	for (i = 0; i=n; i++){
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

	double *At = (double*)calloc(n*n,sizeof(double*));
	double *Y = (double*)calloc(n,sizeof(double));
	//At = D + T*(A-D);
	cblas_dcopy(n,X,1,Y,1);
	APP = PE;
	sc = COUNT(At,PE);
	
	if (sc!=k && sc!=k-1){ReduceStep(A, D, XT, DDD, OO, H, PT, Z, Z1, APP); return 0; }	
	if (sc==k){KK=1;}
	if (sc==k-1){KK=0;}
	
	for (i = 0; i<=10; i++){
		Y = II(At,Y,APP,1);
		//Compute Residual
		if (RES<=EPS1){
			sc = COUNT(At,APP-EPS1-EPS3);
			if (sc==k-1){
				Predict(A, D, X, DDD, OO, H, T, Z, Z1, APP);
				return 0;			
			}
			else if (sc>k-1){
				ReduceStep(A, D, XT, DDD, OO, H, PT, Z, Z1, APP);
				return 0;
			}
			else{
				sc = COUNT(At,APP+EPS1+EPS3);
				if (sc>=k){
					Predict(A, D, X, DDD, OO, H, T, Z, Z1, APP);
					return 0;				
				}
			}		
		}
		
		APP = RQI(At,Y,1);
		
		if (KK==1 && APP>(PE+EPS2)){
			ReduceStep(A, D, XT, DDD, OO, H, PT, Z, Z1, APP);
			return 0;
		}
		
		if (KK==0 && APP<(PE-EPS2)){
			ReduceStep(A, D, XT, DDD, OO, H, PT, Z, Z1, APP);
			return 0;
		}
		
		if (i==10){
			ReduceStep(A, D, XT, DDD, OO, H, PT, Z, Z1, APP);
			return 0;
		}
	}
	return 0;		
}

int ReduceStep(double* A, double* D, double* XT, double* DDD, double* OO double H, double PT, double Z, double Z1, double APP){
	double *X = (double*)calloc(n,sizeof(double));
	double T;
	H = H/2;
	T = PT+H;
	cblas_dcopy(n,XT,1,X,1);
	if (NS==0){
		PE = Z+Z1*H;
		mainblock(A, D, DD, DDD, OO, X, XT, Z, Z1, H, T, PT, PE, NS);
	}
	else{
		Predict(A, D, X, DDD, OO, H, T, Z, Z1, APP);
	}		
	return 0;	
}

int Predict(double* A, double* D, double* X, double* DDD, double* OO, double H, double T, double Z, double Z1, double APP){
	double *XT = (double*)calloc(n,sizeof(double));
	double *HA = (double*)calloc(n,sizeof(double));	
	double OZ=Z, OZ1=Z1, Q, QQ, PH=H, PT=T,PE1;
	H = 1-PT; T = 1; NS+=1; 
	cblas_dcopy(n,X,1,XT,1);
	if (T==1){
		//Store the eigenpair	
	}
	else{
		//Recompute Z1
		//
		//
		Q = pow(1+(PH/H),2);
		QQ = Q*(H/PH);
		PE1=OZ+OZ1*(H+PH)+Q*((Z-OZ)-OZ1*PH)+QQ*(PH*(Z+OZ1)-2*(Z-OZ));
		PE = PE1;
		mainblock(A, D, DD, DDD, OO, X, XT, Z, Z1, H, T, PT, PE, NS);
	}
	return 0
}

/*Function generates the type of matrix as described in Example 1*/
int Example1(double* A, int n, int m, int k){
	int i,j,r;
	double l;
	for (i = 0; i<=n; i++){
		A[i*n+i]= (double) i;
		for (j = 1; j<=n; j++){
			r = j%2;
			if (r==0){
				l = sqrt(j-(0.5*j));
				A[j*n+(j-1)]=l
				A[(j-1)+n*j]=l;
			}
			else{
			    l = k*sqrt(m+(j-1)-(0.5*(j-1)));
				A[j*n+(j-1)]=l;
				A[(j-1)+n*j]=l;
			}						
		}
	}
	return 0;
}

/*Function generates the type of matrix as described in Example 3*/
int Example3(double* A, int n){
	int i,j;
	for(i = 0; i<=n; i++){
		A[i*n+i]= (double) i;
		for (j = 1; j<=n; j++){
			A[j*n+(j-1)]=1.0;
			A[(j-1)+j*n]=1.0;		
		}	
	}
	return 0;
}

/*Function generates an identity matrix of size n*/
int Identity(double* I, int n){
	int i;
	for (i = 0; i<=n; i++){
		I[i*n+i] = 1.0;	
	}
	return 0;
}


/*Initialize all of the vectors needed for computation*/
int Initialize(double *A, double *D, double* DD, double* DDD, double* d, double* OO, int n){
	int i;
	for (i = 0; i<=n; i++){
		DD[i] = A[i*n+i];
		D[i*n+i] = DD[i];
		d[i] = D[i*n+i];
		DDD[i] = DD[i]-d[i];
	}
	OO[0]=0;
	for (i = 1; i<n; i++){
		OO[i]=A[i*n+(i-1)];
	}
	return 0;
}

/*Compute the number of sign changes using Sturm Sequence method*/
int COUNT(double *At, double* DD, double* OO, double x){
	int Count = 0, d = 1,i;
	
	for (i = 0, i<=n, i++){
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
	double *W = (double*)calloc(n,sizeof(double));
	double *I = (double*)calloc(n*n,sizeof(double));
	double *Y = (double*)calloc(n,sizeof(double));
	Identity(I,n);	
	for (int i=0; i<=j; i++){
		//Perform Inverse Iteration	
	}
	free(W);
	free(I);
	return Y;
}

/*Approximate an eigenvalue with Rayleigh Quotient*/
double RQI(double *At, double *X, int j){
	double *u = (double*)calloc(n,sizeof(double));
	double *W = (double*)calloc(n,sizeof(double));
	for (int i = 0; i<=j; i++){
		//Peform Rayleigh Quotient Iteration
	}
	free(u);
	free(w);
	return APP;
}


