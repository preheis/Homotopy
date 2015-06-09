function [Count] = COUNT(At,x)
%COUNT finds the number of eigenvalues less than the predicted eigenvale.
%This ensures that we are on the correct eigenpath.
%At - the matrix at t = T
%x - the predicted eigenvalue

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

