function []= newHomotopy(A,a)
%This program computes all of the eigenvalues of a symmetric tridiagonal
%matrix with well-separated eigenvalues.
%A - the matrix whose eigenvalues we are looking for.

%%%%%%%%%%%%
%%%BLOCK1%%%
%%%%%%%%%%%%

%This block initializes the problem. It determines the kth eigenpair 
%of the initial diagonal matrix D. It also determines the predicted 
%eigenvalue at the next step by using Taylor method. 

global k n OO DDD DD EV EVT

[~,n]=size(A);
EV=zeros(n,1);
EVT=zeros(n,n);
%for j=1:n  
k=a;                                    %The kth eigenpair
D=zeros(n,n);
DD=diag(A);                          %Diagonal of A
OO=zeros(n,1);                       %Off-diagonal of A

  for i=1:n
      D(i,i)=DD(i); 
  end

%  for i=1:n
%     if i==1
%          D(i,i)=(1+(4/n));
%      else
%          D(i,i)=(2+(i*(4/n)));
%      end
%  end

OO(1)=0;

for i=2:n
   OO(i)=A(i-1,i);
end

d=diag(D);                           %Diagonal of D

DDD = zeros(n,1);
for i=1:n
    DDD(i) = DD(i)-d(i);
end

alpha=0;
Z=d(k);
I=eye(n);

XT=zeros(n,1);
XT(k)=1;

% %To predict the eigenvalue we use a third order taylor method. 
 
Z1=XT'*(A-D)*XT;                     %Compute first derivative
XT1=(pinv(Z*I-A)*(A-D))*XT;
Z2=-2*(XT1')*XT1-((XT1')*XT1);       %Compute second derivative
XT2=2*(A-D)*XT1+(alpha*Z1*XT1);
Z3=-3*(XT1')*XT2-(3*(XT1')*XT2);     %Compute third derivative 
 
NS=0;H=1;PT=0;T=1;X=XT;                             
 
PE=Z+(H*Z1)+(H^2/2)*Z2+(H^3/6)*Z3;   %Predict the eigenvalue 
fprintf('The predicted eigenvalue is: %2.5e\n',PE);
mainblock(A,D,Z,Z1,H,T,PT,X,XT,PE,NS);  %Begin correction
%end
end

function[] = taylorEstimation(A,D,Z,X,XT,H,T,PT,NS)
%To predict the eigenvalue we use a third order taylor method. 
global n 
I = eye(n);
alpha = -1;
%At = D+T*(A-D);

disp('Using Taylor Method to predict eigenvale');
Z1=XT'*(A-D)*XT;                     %Compute first derivative
XT1=(pinv(Z*I-A)*(A-D))*XT;
Z2=-2*(XT1')*XT1-((XT1')*XT1);       %Compute second derivative
XT2=2*(A-D)*XT1+(alpha*Z1*XT1);
Z3=-3*(XT1')*XT2-(3*(XT1')*XT2);     %Compute third derivative 
PE=Z+(H*Z1)+(H^2/2)*Z2+(H^3/6)*Z3;   %Predict the eigenvalue 

fprintf('The new predicted eigenvalue is %2.20f\n',PE);
mainblock(A,D,Z,Z1,H,T,PT,X,XT,PE,NS);
end


function [] = mainblock(A,D,Z,Z1,H,T,PT,X,XT,PE,NS)
global k n DD OO
%Set the error tolerances
SUMA=0;
SUMB=0;
DD=diag(A);
for i = 1:n
    SUMA=SUMA+abs(DD(i));
end

for i=1:n
    SUMB=SUMB+abs(OO(i));
end

%NORM is defined as the sum of the sum of diagonal and off diagonal elements 
NORM = SUMA + SUMB;
EPS=(eps/2)*n*NORM;

At = D + T*(A-D);

if T==1
    EPS1=EPS;
else
    dt=1/norm(pinv(Z*eye(n)-At));
    arr=[dt*sqrt(EPS),EPS];
    EPS1=max(arr);
end
EPS2=EPS; EPS3=EPS;
fprintf('The roundoff error is: %2.5e\n',EPS1);

%%%%%%%%%%%%
%%%BLOCK2%%%
%%%%%%%%%%%%

APP=PE;
fprintf('The value of T is: %2.20f\n',T);
sc=Count(At,PE);                               %Compute the number of sign changes
fprintf('The number of sign changes is %d\n',sc);
pause
if (sc~=k)&&(sc~=k-1)
    disp('Reducing the step size!');
    reduceStep(A,D,H,PT,XT,NS,Z,Z1,APP);       %Reduce the step size
    return
end

if sc==k 
    KK=1;
end

if sc==k-1
    KK=0;
end

Y=X;
fprintf('The value of KK is %d\n',KK);
fprintf('The approximate eigenvalue is %2.10f\n',APP);
disp('The predicted eigenvector is: ');
%disp(Y');
pause
for i = 1:10
[Y,RES] = II(At,Y,APP,1);                            %Perform inverse iteration
fprintf('The Residual value is: %2.20e\n',RES);
fprintf('The roundoff error is: %2.20e\n',EPS1);
pause
if RES<=EPS1
   sc=Count(At,APP-EPS1-EPS3);                       %Compute the number of sign changes
   fprintf('The number of sign changes is %d\n',sc);
   if sc==k-1
      predict(A,D,NS,Z,Z1,APP,Y,H,T);                %Store/Recompute
      return
   elseif sc>k-1
      disp('Reducing the step size!');
      reduceStep(A,D,H,PT,XT,NS,Z,Z1,APP);            %Reduce the step size
      return
   else
      sc=Count(At,APP+EPS1+EPS3);
      if sc>=k
         predict(A,D,NS,Z,Z1,APP,Y,H,T);             %Store/Recompute
         return
      end
   end
end

APP=RQI(At,Y,1);                                      %Perform Rayleigh Quotient Iteration
fprintf('The corrected eigenvalue is %2.5e\n',APP);
pause
if KK==1&&(APP>(PE+EPS2))%Check if APP is reasonable
    fprintf('APP is: %2.5e\n',APP);
    fprintf('PE + EPS2 is: %2.5e\n',(PE+EPS2))
    disp('Reducing the step size!');
    reduceStep(A,D,H,PT,XT,NS,Z,Z1,APP);              %Reduce the step size
    return
end

if KK==0&&(APP<(PE-EPS2))                             %Check if APP is reasonable
    fprintf('APP is: %2.5e\n',APP);
    fprintf('PE - EPS2 is: %2.5e\n',(PE-EPS2));
    disp('Reducing the step size!');   
    reduceStep(A,D,H,PT,XT,NS,Z,Z1,APP);              %Reduce the step size
    return
end

if i==10
    disp('Too many iterations of RQ/II');
    reduceStep(A,D,H,PT,XT,NS,Z,Z1,APP);              %Reduce the step size
    return
end
end
end

function[]=reduceStep(A,D,H,PT,XT,NS,Z,Z1,APP)
H=H/2;                                               %Cut step size in half
T=PT+H;                                              %Decrease the homotopy parameter
X=XT;
fprintf('The value of H is %2.5e\n',H);
fprintf('The value of T is %2.5e\n',T);
fprintf('The number of steps is %d\n',NS)
pause
if NS==0
    disp('Generating new prediction using taylor estimation');
    taylorEstimation(A,D,Z,X,XT,H,T,PT,NS);            %Use third order taylor method to predict again
else
    disp('Computing eigenvalue derivative and generating new prediction');
    predict(A,D,NS,Z,Z1,APP,X,H,T);                  %Predict eigenvalue with hermite interpolation or store the result
end
end

function [] = predict(A,D,NS,Z,Z1,APP,X,H,T)
%%%%%%%%%%%%
%%%BLOCK3%%%
%%%%%%%%%%%%

global n k OO DDD
NS=NS+1;OZ=Z;OZ1=Z1;Z=APP;XT=X;
fprintf('The value of T is %2.5e\n',T);
if T==1
    disp('Storing the eigenvalues and eigenvectors!');
    store(Z,XT);
    return
else
    HA=zeros(n,1);
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
fprintf('Z is: %2.5e\n',Z);
fprintf('Z1 is: %2.5e\n',Z1);
fprintf('OZ is: %2.5e\n',OZ);
fprintf('OZ1 is: %2.5e\n',OZ1);
%disp('XT is: ');
%disp(XT);
PH=H;PT=T;H=1-PT;T=1;
fprintf('PT is: %2.5e\n',PT);
fprintf('PH is: %2.5e\n',PH);
fprintf('H is: %2.5e\n',H);
%Use hermite interpolation to get a new prediction.
disp('Using hermite interpolation to predict the eigenvalue!');
Q=(1+(H/PH))^2;
fprintf('Q is: %2.5e\n',Q);
QQ=Q*(H/PH);
fprintf('QQ is: %2.5e\n',QQ);
PE=OZ+OZ1*(H+PH)+Q*((Z-OZ)-OZ1*PH)+QQ*(PH*(Z1+OZ1)-2*(Z-OZ));

fprintf('The value of T is: %2.5e\n',T);
fprintf('The new predicted eigenvalue is: %2.5e\n',PE);
pause
mainblock(A,D,Z,Z1,H,T,PT,X,XT,PE,NS);
end

function [] = store(Z,XT)
global n k EV EVT

fprintf('\n\n')

EV(k)=Z;
for i=1:n
   EVT(i,k)=XT(i);
end

disp(EV);
%disp(EVT);

if k==n
    A=matGen(n,1);
    disp('The eigenvalues are: ')
    disp(EV);
    disp('The eigenvectors are: ')
    %disp(EVT);
    ORT=EVT'*EVT-eye(n);
    disp('The maximum orthogonality is: ');
    disp(max(max(ORT)));
    disp('The maximum residual is: ')
    for i=1:n
        RES = A*EVT(:,i)-EV(i)*EVT(:,i);
    end
    RES = max(RES);
    disp(RES);
end
end

function [Count] = Count(At,x)
%COUNT finds the number of eigenvalues less than the predicted eigenvale.
%This ensures that we are on the correct eigenpath.
%At - the matrix at t = T
%x - the predicted eigenvalue
Count = 0;
d=1;
[~,n]=size(At);
a=diag(At);
b=zeros(n,1);
b(1)=0;
for i=2:n
    b(i)=At(i,i-1);
end

for i=1:n
    d=a(i)-x-(b(i)^2)/d;
    if d<0
        Count=Count+1;
    end
end

end

function [Y,RES]=II(A,X,APP,k)
%Perform inverse iteration to obtain a better approximate eigenvector.
W=A-APP*eye(size(A));
for j=1:k
    Y=X/norm(X);
    X=W\Y;
end
Y=X/norm(X);
RES=A*Y-APP*Y;
RES=abs(max(RES));
end

function [APP]=RQI(A,X,k)
%Perform Rayleigh quotient iteration to obtain a better approximate
%eigenvale. 
for j=1:k
    u=X/norm(X);
    APP=u'*A*u;
    As=A-APP*eye(size(A));
    %fprintf('The condition of the matrix is: %2.20f\n',rcond(As));
    X=(As)\u;
end
APP=u'*A*u;
end