function [] = Homotopy(A,a)
%This function computes the eigenvalues and eigenvectors of a symmetric
%tridiagonal matrix A. Specify the eigenpair a. 

%INITIALIZE HOMOTOPY ALGORITHM
[n,n] = size(A);

global NS T

%D is the diagonal of the intial matrix B
D = zeros(1,n);
%DD and OO are diagonal and off diagonal elements of A
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


%Set the error tolerances
SUMA = 0;
SUMB = 0;

for i = 1:n
    SUMA=SUMA+abs(DD(i));
    SUMB=SUMB+abs(OO(i));
end

%NORM is defined as the sum of the sum of diagonal and off diagonal elements 
NORM = SUMA + SUMB;

EPS = eps*n*NORM;

global EPS1 EPS2 EPS3

if T == 1
    EPS1 = EPS;
else
    %arr = [d0*sqrt(EPS),EPS];
    %EPS1 = max(arr);
    EPS1 = EPS;
end

EPS2 = EPS; EPS3 = EPS;

%__________________100___________________%
function PE = initialPrediction(H)
%Use third order taylor method to predict eigenvalue
    PE = Z + (H*Z1) + (H^2/2)*Z2 + (H^3/6)*Z3;
end

PE = initialPrediction(H);

%__________________200___________________%
function [KK,APP,F,G,J] = main(PE)
%Locate the kth eigenpair at t = T, starting with the predicted eigenvalue
%and using Rayleigh quotient and inverse iteration. If this fails we cut
%the step size in half and start over from when the problem is initialized.

%Form the one parameter family of matricies;
t = T;

At = B + t(A-B);

%J = iterations of RQI
J = 0;
APP = PE;
%F = zeros(1,n);
G = zeros(1,n);

F = D - T*DDD;

for j = 1:n
    F(j) = D(j)*DDD(j);
end

for j = 1:n
    G = T*OO(j);
end

%Computer Sturm Sequence of At 
SC = COUNT(At,PE);

%k is the kth eigenvalue of A(t)
if SC ~= k && SC ~= k-1 
    reduceStep(H,PT);
end

if SC == k
    KK = 1;
end

if SC == k-1
    KK = 0;
end

end

KK = main(PE);

function invit(At,APP,J)
%Here we calculate inverse iteration and normalize the resulting vector.

%W is real symmetric and tridiagonal. 
%RES is the norm of the residual vector @ current eigenpair approximation

I = eye(n);
W = At - APP*I;

%B is not used anywhere else in the algorithm
B = F - APP;

Y = zeros(1,n);

for l = 1:n
    Y(l) = X(l);
end

X = Y\W;

SUM = 0;

for l = 1:n
    SUM = SUM + X(l)^2;
end

RES = 1/sqrt(SUM);

for l = 1:n
    X(l) = X(l)*RES;
end

%Stop correction if residual vector is less than prescribed accuracy
%Else apply Rayleigh Quotient 
if RES<=EPS1
    checkError(APP);
else
    %PAPP isn't used anywhere else in the algorithm. 
    %PAPP = APP;
    
    H(1) = F(1)*X(1)+G(2)*X*(2);
    for l = 2:n-1
        H(l) = G(l)*X(l-1)+F(l)*X(l)+G(l+1)*X(l+1);
    end
    
    H(n) = G(n)*X(n-1)+F(n)*X(n);
    
    APP = 0;
    
    for l = 1:n
        APP = APP+X(l)*H(l);
    end
    
    if (KK == 1) && (APP>(PE+EPS2))
        reduceStep(H,T);
    end
    
    if (KK == 0) && (APP<(PE-EPS2))
        reduceStep(H,T);
    end
    
    J = J+1;
    
    if J == 10
        reduceStep(H,T);
    else
        %recursively call invit 
        invit(At,APP,J);
    end
end

end

%__________________400___________________%
function checkError(APP)

err1 = APP-EPS1-EPS3;

val = COUNT(At,err1);

if val == k-1
    computeDerivative()%GOTO 600
elseif val > k-1
    reduceStep()%GOTO 500  
else
   err2 = APP+EPS1+EPS3;
   val = COUNT(err2);
   if val >= k
       computeDerivative()%GOTO 600
   end
end

end

%__________________500___________________%
function reduceStep(H,PT,NS)
%If predicted eigenvale is rejected, we reduce stepsize, alter parameter
%and redo the prediction

H = H/2; T = PT+H;

for l = 1:n
    X(l) = XT(l);
end

if NS == 0
    initialPrediction(H);%GOTO 100
else
    newPrediction(H,T);%GOTO 700n
end

end

%__________________600___________________%
    function [NS,OZ,OZ1,Z1] = computeDerivative(NS,Z,Z1,APP)
%If the eigepair of At was successfully located in Block 2, we store it if 
%T = 1, otherwise we compute the eigenvalue derivative and predict the
%eigen pair again. 

NS = NS+1; OZ = Z; OZ1 = Z1; Z = APP; 

for m = 1:n
    XT(m) = X(m);
end

HA = zeros(1,n);

%If T = 1 and the eigenvalue was acceptable, we store it.
%We also store the corresponding eigenvector
if T == 1
   storeFinal(Z,XT);
else
    %Compute the derivative
    HA(1) = DDD(1)*X(1)+OO(2)*X(2);
    for m = 2:n-1 
        HA(m) = OO(m)*X(m-1)+DDD(m)*X(m)+OO(m+1)*X(m+1);
    end
    HA(n) = OO(n)*X(n-1)+DDD(n)*X(n);
    
    Z1 = 0;
    
    for m = 1:n
        Z1 = Z1 + X(m)*HA(m);
    end
end

end
%__________________700___________________%
    function [PH,PT,H,T,PE] = newPrediction(H,T)
%H step size
%T homotopy parameter

%If the first prediction fails, we use hermite interpolation, which
%interpolates current and previous eigenvalues together with their
%derivatives. 
PH = H; PT = T; H = 1-PT; T = 1;

Q = (1+(H/PH))^2;
QQ = Q*(H/PH);
PE = OZ + OZ1*(H+PH)+Q*((Z-OZ)-OZ1*PH)+QQ*(PH*(Z1+OZ1)-2*(Z-OZ));

end

%__________________800___________________%
function storeFinal(Z,XT)
%Store the final accepted Eigenvalue and Eigenvector
%Z the accepted eigenvalue
%XT the accept eigenvector
 
EV = zeros(1,n);
EVT = zeros(n,n);
EV(a) = Z;

for p = 1:n
   EVT(p) = XT(p);
end
    
fprintf('%3.10f\n',EV);
fprintf('%3.10f',EVT);

end
%_________________END____________________%
end

%Calculate the number of sign changes using sturm sequence 
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

