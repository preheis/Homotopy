A = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];

a=4;
d = diag(A);
[~,n] = size(A);
D = zeros(n,n);

for i = 1:n
    D(i,i) = d(i);
end

alpha = 1;
%Predict the 1st eigenvalue
Z = d(a);
I = eye(n);

XT = zeros(n,1);

XT(a) = d(1);

Z1 = XT'*(A-D)*XT;
XT1 = pinv(Z*I-A)*(A-D)*XT;
Z2 = -2*(XT1')*XT1-((XT1')*XT1);
XT2 = 2*(A-D)*XT1+(alpha*Z1*XT1);
Z3 = -3*(XT1')*XT2-(3*(XT1')*XT2);
H = 1;

PE = Z+(H*Z1)+(H^2/2)*Z2+(H^3/6)*Z3;

disp(PE)
