function [PE] = RQI(A,x,k)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
for j = 1:k
    u = x/norm(x);
    PE = u'*A*u;
    As = A-PE*eye(size(A));
    if rcond(As)< eps
        break
    end
    x = (As)\u;
end
PE = u'*A*u;
end

