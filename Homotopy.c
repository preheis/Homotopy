#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/*This program computes all of the eigenpairs of symmetric tridiagonal matricies using a Homotopy algorithm*/

//Global variables
int N = 4;
int k = 0;

/*This function reduces the step size, alters the homotopy parameter, and resets the eigenvector*/
int ReduceStep(double *A, double *D, double *DD, double *XT, double *DDD, double *OO, double H, double PT, double Z, double Z1, double APP, int NS){
	double *X = (double*)calloc(N,sizeof(double));
	double T;
	H = H/2;
	T = PT+H;
	cblas_dcopy(N,XT,1,X,1);
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
	double *XT = (double*)calloc(N,sizeof(double));
	double *HA = (double*)calloc(N,sizeof(double));
	double *EV = (double*)calloc(N,sizeof(double));
	double *EVT= (double*)calloc(N*N,sizeof(double));	
	double OZ=Z, OZ1=Z1, Q, QQ, PH=H, PT=T,PE;
	int i;
	H = 1-PT; T = 1; NS+=1; 
	cblas_dcopy(N,X,1,XT,1);
	if (T==1){
		EV[k]=Z;
			
	}
	else{
		//Recompute Z1
		HA[0]=DDD[0]*X[0]+OO[1]*X[1];
		for (i = 1; i<N-1; i++){
			HA[i]=OO[i]*X[i-1]+DDD[i]*X[i]+OO[i+1]*X[i+1];
		}
		HA[N-1]=OO[N-1]*X[N-2]+DDD[N-1]*X[N-1];
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
	for(i = 0; i<N; i++){
		A[i*N+i]= (double) (i+1);
		for (j = 0; j<N; j++){
			A[(j*N)+(j-1)]=1.0;
			A[((j-1)*N)+j]=1.0;
		}	
	}
	
	return 0;
}

/*Function generates an identity matrix of size n*/
int IdentityMat(double *ID){
	int i;
	double j=1;
	for (i = 0; i<N; i++){
		ID[(i*N)+i] = j;	
	}
	return 0;
}


/*Initialize all of the vectors needed for computation*/
int Initialize(double *A, double *D, double *DD, double *DDD, double *OO){
	int i;
	for (i = 0; i<N; i++){
		DD[i] = A[i*N+i];
		D[i*N+i] = DD[i];
		DDD[i] = DD[i]-DD[i];
	}
	OO[0]=0;
	for (i = 1; i<N; i++){
		OO[i]=A[i*N+(i-1)];
	}
	return 0;
}

/*Compute the number of sign changes using Sturm Sequence method*/
int COUNT(double *DD, double *OO, double PE){
	int Count = 0,i;
	double d;
	d = 1.0;
	for (i = 1; i<N; i++){
		d = DD[i]-PE-((OO[i-1]*OO[i-1])/d);
		if (d<0){Count++;}
	}
	return Count;
}

//////////////////////////////////////////////////////////////////////

/*This function computes the residual of the vector obtain from II*/
double computeRES(double* A, double* Y, double APP){
	double *AY = (double*)calloc(N,sizeof(double));
	double *APPY = (double*)calloc(N,sizeof(double));
	double *RES = (double*)calloc(N,sizeof(double));	
	double maxRES=0.0;
	int incx=1;
	int incy=1;
	int alpha=1;
	int beta=1;
	int i;

	cblas_dcopy(N,Y,incy,APPY,incy);//APPY <- Y
	//printMat(A);
	cblas_dscal(N,-APP,APPY,incy);//APPY <- APP*Y
	dgemv_("n",&N,&N,&alpha,A,&N,Y,&incx,&beta,AY,&incy);//AY <- A*Y
	for (i = 0; i<N; i++){RES[i] = AY[i]-APPY[i];}
	printf("Finding the maximal element...\n");
	printVec(RES);
	//maxRES = isamax_(&N,RES,&incx);//RES <- max(RES);
	for (i = 0; i<N; i++){
		if (RES[i]>maxRES){
			maxRES = RES[i];
		}	
	}
	printf("The residual is %.5e\n",maxRES);
	free(AY);
	free(APPY);
	free(RES);
	return maxRES;
}

/*Approximate an eigenvector with Inverse Iteration*/
double * II(double *At, double *X, double APP, int j){
	printf("Beginning Inverse Iteration...\n");
	double *W = (double*)calloc(N*N,sizeof(double));
	double *ID = (double*)calloc(N*N,sizeof(double));
	double *Y = (double*)calloc(N,sizeof(double));
	int *IPIV = (int*)malloc(N*sizeof(int));
	double norm=0, SUM = 0;	
	int i,k,NRHS=1,LDB=N,LDA=N, INFO;
	IdentityMat(ID);
	//W = At-APP*I
	for (i = 0; i<N; i++){
		for (k = 0; k<N; k++){
			W[i*N+k] = At[i*N+k] - APP*ID[i*N+k];		
		}	
	}
	//for (i=0; i<=j; i++){
		//Perform Inverse Iteration
		//norm = cblas_dnrm2(n,X,1); //norm(X);

		for (i=0; i<N; i++){
			SUM += X[i]*X[i];
		}

		norm = sqrt(SUM);

		cblas_dscal(N,1/norm,X,1); // X <- X*1/norm(X)
		printf("Computing the LU factorization...\n");
		dgetrf_(&N,&N,W,&LDA,IPIV,&INFO);
		printf("Solving the system of linear equations...\n");
		dgetrs_("n",&N,&NRHS,W,&LDA,IPIV,X,&LDB,&INFO);			
	//}
	//norm = cblas_dnrm2(n,X,1); //norm(X);
	for (i=0; i<N; i++){SUM += X[i]*X[i];}
	norm = sqrt(SUM);
	cblas_dscal(N,1/norm,X,1); // X <- X*1/norm(X)
	cblas_dcopy(N,X,1,Y,1); // Y <- X
	printVec(Y);	
	free(W);
	free(ID);
	free(IPIV);
	printf("Ending Inverse Iteration\n");
	return Y;
}

/*Approximate an eigenvalue with Rayleigh Quotient*/
double RQI(double *At, double *Y, int j){
	double *u = (double*)calloc(N,sizeof(double));
	double norm,APP,SUM;
	int i;
	int incx=1;
	double alpha=1.0;
	double beta=1.0;
	int incy=1;
	printf("Beginning Rayleigh Quotient Iteration...\n");
	printVec(Y);
	//Peform Rayleigh Quotient Iteration
	//norm = cblas_dnrm2(n,X,incx); //norm(X)
	for (i=0; i<N; i++){SUM += Y[i]*Y[i];}
	norm = sqrt(SUM);
	printf("Norm = %.5f\n",norm);
	cblas_dscal(N,1/norm,Y,incx); //X <- X*1/norm(X);
	printVec(Y);
	printMat(At);
  //dgemv_("n",&N,&N,&alpha,AD,&N,XT,&incx,&beta,XTprime,&incy); 
	dgemv_("n",&N,&N,&alpha,At,&N,Y,&incx,&beta,u,&incy); //u <- At*X
	printVec(u);
	APP = ddot_(&N,Y,&incx,u,&incy); //APP <- X'*u
	printf("Ending Rayleigh Quotient Iteration...\n");
	printf("APP = %.5f\n",APP);
	free(u);
	return APP;
}
/*The main block of the program. It corrects the prediction computed.*/
int mainblock(double *A, double *D, double *DD, double *DDD, double *OO, double *X, double *XT, double Z, double Z1, double H, double T, double PT, double PE, int NS){

	double eps, EPS, EPS1, EPS2, EPS3, NORM, SUMA, SUMB, APP, RES;
	int i,sc,KK,j;
	
	//Set the error tolerances for the program
	eps = 1.11E-16;

	for (i = 0; i<N; i++){
		SUMA+=abs(DD[i]);
		SUMB+=abs(OO[i]);
	}
	
	NORM = SUMA+SUMB;
	
	EPS=eps*N*NORM;
	EPS2=EPS3=EPS;
	if (T==1){
		EPS1=EPS;
	}
	else{
		EPS1=sqrt(EPS);
	}
	
	//////////////////////////////////////////////////

	double *At = (double*)calloc(N*N,sizeof(double));
	double *Y = (double*)calloc(N,sizeof(double));
	
	//At = D + T*(A-D);
	
	for (i=0; i<N; i++){
		for (j=0; j<N; j++){
			At[i*N+j] = D[i*N+j]+T*(A[i*N+j]-D[i*N+j]);		
		}
	}

	cblas_dcopy(N,X,1,Y,1);
	APP = PE;
	sc = COUNT(DD,OO,APP);
	
	printf("The number of sign changes is %d\n",sc);
	
	if (sc!=k && sc!=k-1){ReduceStep(A, D, DD, XT, DDD, OO, H, PT, Z, Z1, APP, NS); return 0; }	
	if (sc==k){KK=1;}
	if (sc==k-1){KK=0;}
	
	for (i = 1; i<=10; i++){
		Y = II(At,Y,APP,1);
		printVec(X);
		printVec(Y);
		//Compute Residual
		printf("Computing the residual...\n");
		RES = computeRES(A,Y,APP);
		printf("The residual is: %e\n",RES);
		printf("Epsilon is: %e\n",EPS1);
		if (RES<=EPS1){
			sc = COUNT(DD,OO,APP-EPS1-EPS3);
			if (sc==k-1){
				Predict(A, D, DD, Y, DDD, OO, H, T, Z, Z1, APP, NS);
				return 0;			
			}
			else if (sc>k-1){
				ReduceStep(A, D, DD, XT, DDD, OO, H, PT, Z, Z1, APP, NS);
				return 0;
			}
			else{
				sc = COUNT(DD,OO,APP+EPS1+EPS3);
				if (sc>=k){
					Predict(A, D, DD, Y, DDD, OO, H, T, Z, Z1, APP, NS);
					return 0;				
				}
			}		
		}
		
		APP = RQI(At,Y,1);
		
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
	return 0;		
}

/*This function prints a matrix*/
int printMat(double *A){
	int i,j;
	for (i=0; i<N; i++){
			for (j=0; j<N; j++){
				printf("%.4f\t",A[i*N+j]);		
			}
			printf("\n");		
		}
	printf("\n");
}

int printVec(double *X){
	int i;
	for (i = 0; i<N; i++){printf("%.4f\n",X[i]);}
	return 0;	
}


//Main
int main(void){
	//Other important variables
	double alpha, beta, H, PT, T, Z, Z1, Z2, Z3, PE;	
	int NS=0,i,incx, incy,nincx;

	PT = 0;
	H = T = 1.0;
	int j;
	incx = incy = 1;
	nincx=-1;
	alpha = beta = 1.0;
	//clock_t t; 
	//Matrices
	double *A = (double*)calloc(N*N,sizeof(double));
	double *D = (double*)calloc(N*N,sizeof(double));
	double *AD = (double*)calloc(N*N,sizeof(double));
	//Vectors	
	//double *d = (double*)calloc(n,sizeof(double));
	double *DD = (double*)calloc(N,sizeof(double));
	double *DDD = (double*)calloc(N,sizeof(double));
	double *OO = (double*)calloc(N,sizeof(double));
	double *X = (double*)calloc(N,sizeof(double));	
	double *XT = (double*)calloc(N,sizeof(double));
	double *XTprime = (double*)calloc(N,sizeof(double));

	Example3(A);
	Initialize(A,D,DD,DDD,OO);
	
	for (i = 0; i<N; i++){
		for (j = 0; j<N; j++){
			AD[i*N+j] = A[i*N+j]-D[i*N+j];		
		}
	}
	//t = clock();
	printf("Beginning Computation\n");
	for (i = 0; i<N; i++){
		k=i;
		Z=DD[k];

		printf("Z = %2.5f\n", Z);	

		XT[k]=1;
		
		cblas_dcopy(N,XT,incx,X,incx);//X <- XT
		
		//Compute first derivative
		printf("Computing Z1\n");
	
		dgemv_("n",&N,&N,&alpha,AD,&N,XT,&incx,&beta,XTprime,&incy); //XT' <- (A-D)*XT
		Z1 = cblas_ddot(0,XT,incx,XTprime,incx);//Z1 <- XT'*((A-D)*XT)
		//Predict eigenvalue with first order taylor method
		PE=Z+Z1*H;
		printf("Beginning Correction.\n");
		mainblock(A, D, DD, DDD, OO, X, XT, Z, Z1, H, T, PT, PE, NS);
	}
	//Free all allocated memory
	free(A);
	free(D);
	free(AD);
	free(DD);
	free(DDD);
	free(OO);
	free(X);
	free(XT);
	free(XTprime);
	//t = clock()-t;
	//fprintf("Total execution time: %f seconds\n", (float)t/CLOCKS_PER_SEC);
	return 0;
}

