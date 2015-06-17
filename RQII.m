function [ lam,u ] = RQII(A,x,s,k)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
As = A-s*eye(size(A));

for i = 1:k
    u = x/norm(x);
    x = As\u;
    u=x/norm(x);
    u=x/norm(x);
    
    lam = u'*A*u;
    x = (As)\u;
    u = x*norm(x);
end

fprintf('%2.10f\n',lam);
for i = 1:size(A)
    fprintf('%2.10f\n',u(i));
end

end


