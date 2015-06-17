function [Count] = COUNT(At,x)
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

for i = 1:n
    d = a(i)-x-(b(i)^2)/d;
    if d < 0
        Count = Count+1;
    end
end

end

