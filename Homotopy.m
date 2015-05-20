function [EV,EVT ] = Homotopy(A,a)
%This function computes the eigenvalues and eigenvectors of a symmetric
%tridiagonal matrix A. Specify the eigenpair a. 

%%%%%%%%%%
%%BLOCK1%%
%%%%%%%%%%

%INITIALIZE HOMOTOPY ALGORITHM
[n,n] = size(A);

%D is the diagonal of the intial matrix B
D = zeros(1,n);
%DD and OO are diagonal and off diagonal elements of A
DD = zeros(1,n);
OO = zeros(1,n);

DD = diag(A);
OO = diag(A,-1);
OO(1) = 0;

%B is the initial matrix
B = zeros(n,n);

for i = 1:n
    B(i,i) = DD(i);
end

DDD = zeros(1,n);

for i = 1:n
    DDD(i) = DD(i) - D(i);
end

%Z is the kth smallest eigenvalue and XT is the corresponding eigenvector.
Z = D(a);
XT = zeros(1,n);
XT(a) = 1;

%Z1, Z2, Z3 are the 1st, 2nd, and 3rd derivatives of Z. 
Z1 = XT*(A-B)*XT';
Z2 = 0;
Z3 = 0;

%V1 and V2 are the eigenvector derivatives of 

%NS = number of steps
%H = step size
%T = homotopy parameter
%PE = predicted eigenvalue 

NS = 0; H = 1; PT = 0; T = 1;
X = zeros(1,n);

for i = 1:n
    X(i) = XT(i);
end

%Use third order taylor method to predict eigenvalue
PE = Z + (H*Z1) + (H^2/2)*Z2 + (H^3/6)*Z3;

%%%%%%%%%%
%%BLOCK2%%
%%%%%%%%%%

%Form the one parameter family of matricies;
t = T;

At = B + t(A-B);

%J = iterations of RQI
J = 0;
APP = PE;
F = zeros(1,n);
G = zeros(1,n);

F = D - T*DDD;

for i = 1:n
    G = T*OO(i);
end

%Computer Sturm Sequence of At 
SC = COUNT(At,PE);

%k is the kth eigenvalue of A(t)
if SC ~= k && SC ~= k-1 
    %GOTO 500    
    if NS == 0
        %%GOTO 100
    else
        %GOTO 700
    end
end

if SC == k
    KK = 1;
end

if SC == k-1
    KK = 0;
end


%W is real symmetric and tridiagonal. 
%RES is the norm of the residual vector @ current eigenpair approximation

I = eye(n);
W = At - APP*I;

%B is not used anywhere else in the algorithm
B = F - APP;

Y = zeros(1,n);

for i = 1:n
    Y(i) = X(i);
end

%WX = Y. Doesnt say what to solve for.
product = W*X;
Y = product;

SUM = 0;

for i = 1:n
    SUM = SUM + X(i)^2;
end

RES = 1/sqrt(SUM);

for i = 1:n
    X(i) = X(i)*RES;
end

SUMA = 0;
SUMB = 0;

for i = 1:n
    SUMA=SUMA+abs(DD(i));
    SUMB=SUMB+abs(OO(i));
end

%NORM is defined as the sum of the sum of diagonal and off diagonal elements 
NORM = SUMA + SUMB;

EPS = eps*n*NORM;

if T ==1
    EPS1 = EPS;
else
    %arr = [d0*sqrt(EPS),EPS];
    %EPS1 = max(arr);
end

EPS2 = EPS; EPS3 = EPS;

%Stop correction if residual vector is less than prescribed accuracy
%Else apply Rayleigh Quotient 
if RES<=EPS1
    %GOTO 400
else
    PAPP = APP;
    
    H(1) = F(1)*X(1)+G(2)*X*(2);
    for i = 2:n-1
        H(i) = G(i)*X(i-1)+F(i)*X(i)+G(i+1)*X(i+1);
    end
    H(n) = G(n)*X(n-1)+F(n)*X(n);
    
    APP = 0;
    
    for i = 1:n
        APP = APP+X(i)*H(i);
    end
    
    if (KK == 1) && (APP>(PE+EPS2))
        %GOTO500
    end
    
    if (KK == 0) && (APP<(PE-EPS2))
        %GOTO500
    end
    
    J = J+1;
    
    if J == 10
        %GOTO500
    else
        %GOTO300
    end
    
end

val = COUNT(At,APP-EPS1-EPS3);

if val == k-1
    %GOTO 600
elseif val > k-1
    %GOTO 500    
else
   val = COUNT(APP+EPS1+EPS3);
   if val >= k
       %GOTO 600
   end
end

%If predicted eigenvale is rejected, we reduce stepsize, alter parameter
%and redo the prediction

H = H/2; T = PT+H;

for i = 1:n
    X(i) = XT(i);
end

if NS == 0
    %GOTO 100
else
    %GOTO 700
end


%%%%%%%%%%
%%BLOCK3%%
%%%%%%%%%%

%If the eigepair of At was successfully located in Block 2, we store it if 
%T = 1, otherwise we compute the eigenvalue derivative and predict the
%eigen pair again. 

NS = NS+1; OZ = Z; OZ1 = Z1; Z = APP; 

for i = 1:n
    XT(i) = X(i);
end


%If T = 1 and the eigenvalue was acceptable, we store it.
%We also store the corresponding eigenvector
if T == 1
   
    EV = zeros(1,n);
    EVT = zeros(n,n);
    
    EV(a) = Z;
    
    for i = 1:n
        EVT(i,:) = XT(i);
    end
    
else
    %Compute the derivative
    H(1) = DDD(1)*X(1)+OO(2)*X(2);
    for i = 2:n-1 
        H(i) = OO(i)*X(i-1)+DDD(i)*X(i)+OO(i+1)*X(i+1);
    end
    H(n) = OO(n)*X(n-1)+DDD(n)*X(n);
    
    Z1 = 0;
    
    for i = 1:n
        Z1 = Z1 + X(i)*H(i);
    end
    
end

%If the first prediction fails, we use hermite interpolation, which
%interpolates current and previous eigenvalues together with their
%derivatives. 
PH = H; PT = T; H = 1-PT; T = 1;

Q = (1+(H/PH))^2;
QQ = Q*(H/PH);
PE = OZ + OZ1*(H+PH)+Q*((Z-OZ)-OZ1*PH)+QQ*(PH*(Z1+OZ1)-2*(Z-OZ));
%GOTO 200


end


function SS = COUNT(At,lambda)
    %At the tridiagonal matrix formed at step t
    %lambda the predicted eigenvalue
    %p the polynomial vector
    %s the number of sign changes
    
    n = size(At,1);
    p = zeros(n,1);
    SS = 0;
    
    p(1) = lambda - At(1,1);
    
    if p(1) < 0
        SS = 1;
    end
    
    p(2) = (lambda-T(2,2))*p(1) - abs(T(1,2))^2;
    
    if p(2)*p(1) < 0
        SS = SS+1;
    end
    
    for i = 3:n
        p(i) = (lambda-At(i,i))*p(i-1) - abs(At(i-1,i))^2*p(k-2);
        if p(k)*p(k-1) < 0
            SS = SS+1;
        end
    end
end

