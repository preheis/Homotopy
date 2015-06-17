function [] = testHomotopy(A)
[~,n] = size(A);

for a = 1:n

    d = diag(A);
    D = zeros(n,n);

    for i = 1:n
        D(i,i) = d(i);
    end

    alpha = 1;
    %Predict the 1st eigenvalue
    Z = d(a);
    I = eye(n);

    XT = zeros(n,1);

    XT(a) = d(a);

    Z1 = XT'*(A-D)*XT;
    XT1 = pinv(Z*I-A)*(A-D)*XT;
    Z2 = -2*(XT1')*XT1-((XT1')*XT1);
    XT2 = 2*(A-D)*XT1+(alpha*Z1*XT1);
    Z3 = -3*(XT1')*XT2-(3*(XT1')*XT2);
    H = 1;

    PE = Z+(H*Z1)+(H^2/2)*Z2+(H^3/6)*Z3;

    for i = 1:10
        if i == 1
            x = II2(A,XT,PE,1);
        else
            x = II2(A,x,PE,1);
        end
        PE = RQI(A,x,1);
    end

    fprintf('Z = %2.15f\n\n',PE);
    disp('X =');
    for i = 1:n
        fprintf('%2.15f\n',x(i))
    end
    fprintf('\n')
end

end