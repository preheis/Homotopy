function [EV,EVT] = Homotopy(A,a,E,X)
%This function computes the eigenvalues and eigenvectors of a symmetrical
%tridiagonal matrix.
disp(nargin);
if nargin < 1
    error('Not enough arguments!');
end
global k n OO DD DDD
if nargin == 2
k = a;
[~,n] = size(A);
I = eye(n);
D = zeros(n,n);

%DD and OO are the diagonal and off diagonal elements of A
DD = zeros(n,1);
OO = zeros(n,1);

for i = 1:n
    DD(i) = A(i,i);
end

for i = 1:n
    D(i,i) = DD(i);
end

%d = diag(D);

OO(1) = 0;

for i = 2:n
    OO(i) = A(i,i-1);
end

for i = 1:n
    DDD(i) = DD(i)-OO(i);
end

%Compute the 1st 2nd and 3rd eigenvalue derivatives.
Z = DD(k);
XT = zeros(n,1);
XT(k) = 1;
alpha = 1;

Z1 = XT'*(A-D)*XT;
XT1 = pinv(Z*I-A)*(A-D)*XT;
Z2 = -2*(XT1')*XT1-((XT1')*XT1);
XT2 = 2*(A-D)*XT1+(alpha*Z1*XT1);
Z3 = -3*(XT1')*XT2-(3*(XT1')*XT2);

%NS is the number of steps
%H is the step size
%T is the homotopy paramemter
%PE is the predicted eigenvalue
NS = 0; H = 1; PT = 0; T = 1; X = zeros(n,1);

for i = 1:n
    X(i) = XT(i);
end

PE = Z+(H*Z1)+(H^2/2)*Z2+(H^3)*Z3;
disp('THE PREDICTED EIGENVALUE IS: ')
disp(PE);
disp('CALLING THE MAIN BLOCK');
mainblock(A,D,PE,XT,PT,T,NS,H,Z,Z1);
end

if nargin == 4
    EV = E;
    EVT = X;
    disp(EV);
    disp(EVT);
end

end

function [] = prediction(A,D,Z,PT,T,NS,H,X)
global n
disp('BEGINNING PREDICTION AGAIN');
alpha = 1;
I = eye(n);

XT = X;
Z1 = XT'*(A-D)*XT;
XT1 = pinv(Z*I-A)*(A-D)*XT;
Z2 = -2*(XT1')*XT1-((XT1')*XT1);
XT2 = 2*(A-D)*XT1+(alpha*Z1*XT1);
Z3 = -3*(XT1')*XT2-(3*(XT1')*XT2);

PE = Z+(H*Z1)+(H^2/2)*Z2+(H^3)*Z3;

mainblock(A,D,PE,XT,PT,T,NS,H,Z,Z1);

end


function [] = mainblock(A,D,PE,XT,PT,T,NS,H,Z,Z1)

global n k OO DD DDD 

%_________________________________________________________________________%
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

if T == 1
    EPS1 = EPS;
else
    dt = norm(pinv(eye(n)-D));
    arr = [dt*sqrt(EPS),EPS];
    EPS1 = max(arr);
end

EPS2 = EPS; EPS3 = EPS;
%_________________________________________________________________________%

APP = PE;
X = XT;

F = zeros(n,1);
G = zeros(n,1);

for i = 1:n
    F(i) = D(i,i)+T*DDD(i);
    G(i) = T*OO(i);
end
disp('The value of T is: ');
disp(T);
t = T;

At = D + t*(A-D);

%Calculate the number of changes of the Sturm sequence
SC = Count(At,PE);
disp('THE NUMBER OF SIGN CHANGES IS: ')
disp(SC);

if (SC~=k) && (SC~=k-1)
    disp('BEGIN REDUCING STEP')
    reduceStep(A,D,Z,Z1,APP,NS,H,PT,XT);%GOTO 500
end

if SC == k
    KK = 1;
end

if SC == k-1
   KK = 0;
end

I = eye(n);
W = At-APP*I;
disp(W);

B = zeros(n,1);
Y = zeros(n,1);

for i = 1:n
    B(i) = F(i)-APP;
    Y(i) = X(i);
end
disp(B);
disp(Y);
%Perform inverse iteration with shift.
X = linsolve(W,Y);
SUM = 0;
disp(X);
for i = 1:n
    SUM = SUM+X(i)^2;
end
disp(SUM);
RES = 1/sqrt(SUM);

for i = 1:n
    X(i) = X(i)*RES;
end
disp(X);
pause
%If the norm of the residual vector is smaller than the tolerance,
%correction is stopped. If not, a new approximation for the eigenvalue is
%computed using the Rayleigh Quotient.
disp(RES);
disp(EPS1);
pause
if RES <= EPS1
    %GOTO 400
    disp('NORM OF THE RESIDUAL VECTOR IS SMALLER THAN TOLERANCE')
    val = APP-EPS1-EPS3;
    disp(val);
    sc = Count(At,val);
    disp('THE NUMBER OF SIGN CHANGES IS: ');
    disp(sc);
    if sc == k-1
        disp('BEGINNING COMPUTE AND PREDICT')
        computeANDpredict(A,D,NS,Z,Z1,APP,X,T,H);%GOTO 600
    end
    
    if sc > k-1
        disp('REDUCING STEP SIZE');
        reduceStep(A,D,Z,Z1,APP,NS,H,PT,XT);%GOTO 500
    end
    
    if sc < k-1
        val = APP+EPS1+EPS3;
        disp('COMPUTING SIGN CHANGE AGAIN');
        sc = Count(At,val);
        if sc >= k
            disp('BEGINNING COMPUTE AND PREDICT')
            computeANDpredict(A,D,NS,Z,Z1,APP,X,T,H)%GO TO 600
        end
    end
    
end

%PAPP is never actually used. 
%PAPP = APP;
HA = zeros(n,1);
disp('BEGINNING RAYLEIGH QUOTIENT ITERATION')
pause
for j = 1:10
    
    if j == 10
       disp('TOO MANY APPLICATIONS OF RQI')
       disp('BEGIN REDUCING STEP SIZE')
       reduceStep(A,D,Z,Z1,APP,NS,H,PT,XT);%GOTO 500
    end
    
    HA(1) = F(1)*X(1)+G(2)*X(2);
    for i = 2:n-1
        HA(i) = G(i)*X(i-1)+F(i)*X(i)+G(i+1)*X(i+1);
    end
    HA(n) = G(n)*X(n-1)+F(n)*X(n);

    APP = 0;
    for i = i:n
        APP = APP+(X(i)*HA(i));
    end
    disp('KK is: ');
    disp(KK);
    disp('APP is: ');
    disp(APP);
    disp('The predicted eigenvalue is: ');
    disp(PE);
    disp('EPS2 is: ');
    disp(EPS2);
    pause
    if (KK == 1) & (APP>(PE+EPS2))
        disp('BEGIN REDUCING STEP SIZE')
        reduceStep(A,D,Z,Z1,APP,NS,H,PT,XT);%GOTO 500
    end
    
    if (KK == 0) & (APP<(PE-EPS2))
        disp('BEGIN REDUCING STEP SIZE')
        reduceStep(A,D,Z,Z1,APP,NS,H,PT,XT);%GOTO 500
    end
end

end

function [] = reduceStep(A,D,Z,Z1,APP,NS,H,PT,XT)
%If the predicted or corrected eigenvalue was rejected, we reduce the step
%size and redo the predicition.
global n
    disp('REDUCING STEP SIZE');
    disp('The number of steps is: ')
    disp(NS);
    disp('The step size is: ')
    disp(H);
    disp('PT is: ');
    disp(PT);
    disp('The Vector is: ')
    disp(XT);
    H = H/2;
    T = PT+H;
    X = zeros(n,1);
    disp('The step size is now: ')
    disp(H);
    disp('The value of T is: ')
    disp(T);
    for i = 1:n
        X(i) = XT(i);
    end
    pause
    if(NS==0)
        disp('STEP 0. GOING BACK TO PREDICTION');
        prediction(A,D,Z,PT,T,NS,H,X)%GOTO 100
    else
        disp('STEP > 0. BEGINNING COMPUTE AND PREDICT');
        computeANDpredict(A,D,NS,Z,Z1,APP,X,T,H);
    end
    
end

function [] = computeANDpredict(A,D,NS,Z,Z1,APP,X,T,H)
%Recompute the eigenvalue derivative and then make a new prediction.
global n k OO DDD

NS = NS+1; OZ = Z; OZ1 = Z1; Z = APP;
XT = zeros(n,1);
HA = zeros(n,1);

for i = 1:n
    XT(i) = X(i);
end
disp('The eigenvector is: ');
disp(XT);
disp('THE VALUE OF T IS:')
disp(T);
pause
%If t == 1, we store the eigenpair
if T == 1
    [EV,EVT] = store(Z,XT);
    Homotopy(A,k,EV,EVT);
end
%If t ~= 1, we recompute the eigenvalue derivative.
if T < 1
    HA(1) = DDD(1)*X(1)+OO(2)*X(2);
    for i = 2:n-1
        HA(i) = OO(i)*X(i-1)+DDD(i)*X(i)+OO(i+1)*X(i+1);
    end
    HA(n) = OO(n)*X(n-1)+DDD(n)*X(n);
    Z1 = 0;
    for i = 1:n
        Z1 = Z1+X(i)*HA(i);
    end
end
disp('The new eigenvalue derivative is: ')
disp(Z1);
%Now we use hermite interpolation to obtain a new eigenvalue prediction.
PH = H; PT = T; H = 1-PT; T = 1;
Q = (1+H/PH);
QQ = Q*(H/PH);
PE = OZ+OZ1*(H+PH)+Q*((Z-OZ)-OZ*PH)+QQ*(PH*(Z1+OZ1)-2*(Z-OZ));
disp('RETURNING TO THE MAIN BLOCK');
disp(PE);
mainblock(A,D,PE,XT,PT,T,NS,H,Z,Z1);%GOTO 200

end

%if t = 1, we store the kth eigenpair of A 
function [EV,EVT] = store(Z,XT)
global n k
disp('BEGIN STORING THE EIGENPAIR');
EV = zeros(1,n);
EVT = zeros(n,n);

EV(k) = Z;

for i = 1:n 
    EVT(i) = XT(i); 
end

end

function [Count] = Count(At,x)
%COUNT finds the number of eigenvalues less than the predicted eigenvale.
%This ensures that we are on the correct eigenpath.
%At - the matrix at t = T
%x - the predicted eigenvalue
disp('COMPUTE SIGN CHANGES');
 Count = 0;
 d = 1;
 a = diag(At);
 b = diag(At,-1);
 [~,n] = size(At);
 
 for i = 1:n-1
     if i == 1
         d = a(i)-x;
     else
         d = a(i)-x-(b(i)^2)/d;
     end
     
     if d < 0
         Count = Count+1;
     end
 end

end
