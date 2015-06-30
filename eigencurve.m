function [ ] = eigencurve(A,D)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

[n,n] = size(A);

    for t = 0:.001:1

        At = D + t*(A-D);

        [~,Eig] = eig(At);
        
        Eig = diag(Eig);

        Eig = sort(Eig);
     
        for i = 1:n

            e = Eig(i,1);
            
            hold on;

            plot(t,e);
        end

    end
end


