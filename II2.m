function [X] = II2(A,X,APP,k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
W = A-APP*eye(size(A));

for j = 1:k
    Y = X/norm(X);
    if rcond(W) < eps
        break
    end
    X = W\Y;
end


end

