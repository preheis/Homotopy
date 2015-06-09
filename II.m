
[~,n] = size(A);
DD = diag(A);
OO = zeros(n,1);
for i = 1:n
    if i == 1
        OO(i) = 0;
    else
        OO(i) = A(n,n-1);
    end
end
d = diag(D);

J = 0;
APP = PE;
F = zeros(n,1);
G = zeros(n,1);

for i = 1:n
    DDD(i) = DD(i) - d(i);
end

for i = 1:n
   F(i) = D(i)+(T*DDD(i));
   G(i) = T*OO(i);
end

T = 1;

At = D + T*(A-D);

B = zeros(n,1);
Y = zeros(n,1);
X = XT;

for i = 1:n
    B(i) = F(i)-APP;
    Y(i) = X(i);
end

X = W\Y;

SUM = 0;

for i = 1:n
    SUM = SUM+(X(i)^2);
end

RES = 1/sqrt(SUM);

for i = 1:n
    X(i) = X(i)*RES;
end

HA = zeros(n,1);

PAPP = APP;

HA(1) = F(1)*X(1)+G(2)*X(2);

for i = 2:n-1
    HA(i) = G(i)*X(i-1)+F(i)*X(i)+G(i+1)*X(i+1);
end

HA(n) = G(n)*X(n-1)+F(n)*X(n);

APP = 0;

for i = 1:n
    APP = APP+X(i)*HA(i);
    disp(APP);
end
