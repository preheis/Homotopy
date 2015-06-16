function []= newHomotopy(A,a)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here

%%%%%%%%%%%%
%%%BLOCK1%%%
%%%%%%%%%%%%

%This block initializes the problem. It determines the kth eigenpair 
%of the initial diagonal matrix D. It also determines the predicted 
%eigenvalue at the next step by using Taylor method. 

global k n OO DDD
[~,n]=size(A);
k=a;                                 %The kth eigenpair
D=zeros(n,n);
DD=diag(A);                          %Diagonal of A
OO=diag(A,-1);                       %Off-diagonal of A

for i=1:n
    D(i,i)=DD(i); 
end
d=diag(D);                           %Diagonal of D

DDD = zeros(n,1);
for i=1:n
    DDD(i) = DD(i)-d(i);
end

alpha=1;
Z=d(k);
I=eye(n);

XT=zeros(n,1);
XT(k)=d(1);

%To predict the eigenvalue we use a third order taylor method. 

Z1=XT'*(A-D)*XT;                     %Compute first derivative
XT1=pinv(Z*I-A)*(A-D)*XT;
Z2=-2*(XT1')*XT1-((XT1')*XT1);       %Compute second derivative
XT2=2*(A-D)*XT1+(alpha*Z1*XT1);
Z3=-3*(XT1')*XT2-(3*(XT1')*XT2);     %Compute third derivative 

NS=0;                                %Number of steps
H=1;                                 %Step size
PT=0;                                %Partial T
T=1;                                 %Homotopy paramemter
X=XT;                                %

PE=Z+(H*Z1)+(H^2/2)*Z2+(H^3/6)*Z3;   %Predict the eigenvalue 
mainblock(A,D,Z,Z1,H,T,PT,X,PE,NS);  %Begin correction
end


function [] = mainblock(A,D,Z,Z1,H,T,PT,Y,PE,NS)
global k n 
%Set the error tolerances
SUMA=0;
SUMB=0;
DD=diag(A);
OO=diag(A,-1);
for i = 1:n
    SUMA=SUMA+abs(DD(i));
end

for i=1:n-1
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

%%%%%%%%%%%%
%%%BLOCK2%%%
%%%%%%%%%%%%

APP=PE;
At = D + T*(A-D);
sc=Count(At,PE);                               %Compute the number of sign changes
fprintf('The number of sign changes is %d',sc);

if (sc~=k)&&(sc~=k-1)
    reduceStep(H,PT,XT,NS);                    %Reduce the step size
end

if sc==k 
    KK=1;
end

if sc==k-1
    KK=0;
end
fprintf('APP %2.10f',APP);
fprintf('The number of sign changes is %d',sc);
fprintf('The number of sign changes is %d',sc);

for i = 1:10
[X,RES] = II(A,Y,APP,1);                       %Perform inverse iteration
if RES<=1e-10
   sc=Count(APP-EPS1-EPS3);                    %Compute the number of sign changes
   if sc==k-1
      predict(A,D,NS,Z,Z1,APP,X,H,T);          %Store/Recompute 
   elseif sc>k-1 
      reduceStep(A,D,H,PT,XT,NS,Z,Z1,APP);     %Reduce the step size 
   else
      sc=Count(APP+EPS1+EPS3);
      if sc>=k
         predict(A,D,NS,Z,Z1,APP,X,H,T);       %Store/Recompute
      end
   end
end

PE=RQI(A,X,1);                                 %Perform Rayleigh Quotient Iteration
APP=PE;
Y=X;

if KK==1&&(APP>(PE+EPS2))                      %Check if APP is reasonable
    reduceStep(A,D,H,PT,XT,NS,Z,Z1,APP);       %Reduce the step size
end

if KK==0&&(APP<(PE-EPS2))                      %Check if APP is reasonable 
    reduceStep(A,D,H,PT,XT,NS,Z,Z1,APP);       %Reduce the step size
end

if i==10
    reduceStep(A,D,H,PT,Y,NS,Z,Z1,APP);        %Reduce the step size
end
end
end

function[]=reduceStep(A,D,H,PT,XT,NS,Z,Z1,APP)
H=H/2;                                         %Cut step size in half
T=PT+H;                                        %Increase the homotopy parameter
X=XT;

if NS==0
    %GOTO 100
else
    predict(A,D,NS,Z,Z1,APP,X,H,T);            %Predict eigenvalue with hermite interpolation or store the result
end
end

function [] = predict(A,D,NS,Z,Z1,APP,X,H,T)
%%%%%%%%%%%%
%%%BLOCK3%%%
%%%%%%%%%%%%

global n OO DDD
NS=NS+1;OZ=Z;OZ1=Z1;Z=APP;XT=X;

if T==1
    store(Z,XT);
else
    HA = zeros(n,1);
    HA(1)=DDD(1)*X(1)+OO(2)*X(2);
    for i=2:n-1
        HA(i)=OO(i)*X(i-1)+DDD(i)*X(i)+OO(i+1)*X(i+1);
    end
    HA(n)=OO(n)*X(n-1)+DDD(n)*X(n);
    Z1=0;
    for i=1:n
        Z1=Z1+X(i)*HA(i);
    end
end

PH=H;PT=T;H=1-PT;T=1;
%Use hermite interpolation to get a new prediction.
Q=(1+(H/PH))^2;
QQ=Q*(H/PH);
PE=OZ+OZ1*(H+PH)+Q*((Z-OZ)-(OZ1*PH))+QQ*(PH*(Z1+OZ1)-2*(Z-OZ));

mainblock(A,D,Z,Z1,H,T,PT,X,PE,NS);
end

function [EV,EVT] = store(Z,XT)
global k n
EV=zeros(n,1);
EVT=zeros(n,n);

EV(k)=Z;
for i=1:n
   EVT(i,k)=XT(i);
end
end

function [Count] = Count(At,x)
%COUNT finds the number of eigenvalues less than the predicted eigenvale.
%This ensures that we are on the correct eigenpath.
%At - the matrix at t = T
%x - the predicted eigenvalue
Count=0;
d=1;
a=diag(At);
b=diag(At,-1);
[~,n]=size(At);
 
for i=1:n-1
    if i==1
        d=a(i)-x;
    else
        d=a(i)-x-(b(i)^2)/d;
    end
     
    if d<0
        Count=Count+1;
    end
end
end

function [Y,RES] = II(A,X,APP,k)
%Perform inverse iteration to obtain a better approximate eigenvector.
W=A-APP*eye(size(A));
for j=1:k
    Y=X/norm(X);
    if rcond(W)<eps
        break
    end
    X=W\Y;
end
Y=X/norm(X);
disp(Y);
RES=1/norm(X);
disp(RES);
end

function [PE] = RQI(A,x,k)
%Perform Rayleigh quotient iteration to obtain a better approximate
%eigenvale. 
for j = 1:k
    u = x/norm(x);
    PE = u'*A*u;
    As = A-PE*eye(size(A));
    disp(rcond(As));
    if rcond(As)< eps
        break
    end
    x = (As)\u;
end
PE = u'*A*u;
disp(PE);
end