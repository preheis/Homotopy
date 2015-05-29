function [ output_args ] = Homotopy(A,a)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

global k n OO DDD EPS1 EPS2 EPS3

k = a;

[n,n] = size(A);

%DD and OO are the diagonal and off diagonal elements of A
DD = zeros(1,n);
OO = zeros(1,n);

for i = 1:n
    DD(i) = A(i,i);
end

OO(1) = 0;

for i = 2:n
    OO(i) = A(i,i-1);
end

for i = 1:n
    DDD(i) = DD(i)-OO(i);
end

%Compute the 1st 2nd and 3rd eigenvalue derivatives.
Z = DD(k);
XT = zeros(1,n);
XT(k) = 1;

Z1 = 0;
Z2 = 0;
Z3 = 0;

%NS is the number of steps
%H is the step size
%T is the homotopy paramemter
%PE is the predicted eigenvalue
NS = 0; H = 1; PT = 0; T = 1; X = zeros(1,n);

for i = 1:n
    X(i) = XT(i);
end

PE = Z+(H*Z1)+(H^2/2)*Z2+(H^3)*Z3;

mainblock(PE);

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
    %arr = [d0*sqrt(EPS),EPS];
    %EPS1 = max(arr);
    EPS1 = EPS;
end

EPS2 = EPS; EPS3 = EPS;

end

function [] = mainblock(NS,PE)

APP = PE;
F = zeros(1,n);
G = zeros(1,n);

for i = 1:n
    F(i) = D(i)+T*DDD(i);
    G(i) = T*OO(i);
end

At = D + t*(A-D);

%Calculate the number of changes of the Sturm sequence
SC = Count(At,PE);

if SC~=k && SC~= k-1
    %GOTO 500
end

if SC == k
    KK = 1;
end

if SC == k-1
   KK = 0;
end

I = eye(n);
W = At-APP*I;

B = zeros(1,n);
Y = zeros(1,n);

for i = 1:n
    B(i) = F(i)-APP;
    Y(i) = X(i);
end

%Perform inverse iteration with shift.
X = Y\W;

SUM = 0;

for i = 1:n
    SUM = SUM+X(i)^2;
end

RES = 1/sqrt(SUM);

for i = 1:n
    X(i) = X(i)*RES;
end

%If the norm of the residual vector is smaller than the tolerance,
%correction is stopped. If not, a new approximation for the eigenvalue is
%computed using the Rayleigh Quotient. 

if RES <= EPS1
    %GOTO 400
    val = APP-EPS1-EPS3;
    
    sc = COUNT(At,val);
    
    if sc == k-1
        %GOTO 600
    end
    
    if sc > k-1
        %GOTO 500
    end
    
    if sc < k-1
        val = APP+EPS1+EPS3;
        sc = COUNT(At,val);
        if sc >= k
            %GO TO 600
        end
    end
    
end

PAPP = APP;
HA = zeros(1,n);

for j = 1:10
    
    if j == 10
       reduceStep(NS,H,PT,XT);%GOTO 500
    end
    
    HA(1) = F(1)*X(1)+G(2)*X(2);
    for i = 2:n-1
        HA(i) = G(i)*X(i-1)+F(i)*X(i)+G(i+1)*X(i+1);
    end
    HA(n) = G(n)*X(n-1)+F(n)*X(n);

    APP = 0;
    for i = i:n
        APP = APP+X(i)*HA(i);
    end

    if KK == 1 && (APP>(PE+EPS2))
        reduceStep(NS,H,PT,XT);%GOTO 500
    end

    if KK == 0 && ((APP<PE-EPS2))
        reduceStep(NS,H,PT,XT);%GOTO 500
    end
end

end

%If the predicted or corrected eigenvalue was rejected, we reduce the step
%size and redo the predicition.
function [] = reduceStep(NS,H,PT,XT)
    X = zeros(1,n);
    H = H/2;
    T = PT+H;
    for i = 1:n
        X(i) = XT(i);
    end
    
    if(NS==0)
        %GOTO 100
    else
       computeANDpredict(NS,Z,Z1,APP,H,T,X);
    end
    
end

function [] = computeANDpredict(NS,Z,Z1,APP,H,T,X)

NS = NS+1; OZ = Z; OZ1 = Z1; Z = APP;
HA = zeros(1,n);

for i = 1:n
    XT(i) = X(i);
end

%If t == 1, we store the eigenpair
if T == 1
    store(OZ,XT);
%If t ~= 1, we recompute the eigenvalue derivative.
else
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

%Now we use hermite interpolation to obtain a new eigenvalue prediction.
PH = H; PT = T; H = 1-PT; T = 1;
Q = (1+H/PH);
QQ = Q*(H/PH);
PE = OZ+OZ1*(H+PH)+Q*((Z-OZ)-OZ*PH)+QQ*(PH*(Z1+OZ1)-2*(Z-OZ));

mainblock(PE)%GOTO 200

end

%if t = 1, we store the kth eigenpair of A 
function [] = store(Z,XT)

EV = zeros(1,n);
EVT = zeros(n,n);

EV(k) = Z;

for i = 1:n 
    EVT(i) = XT(i); 
end

fprintf('%3.10f\n',EV);
fprintf('%3.10f',EVT);

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
