function [A] = matGen(n,m,k)
%matGen generates a matrix with given parameters. 
%n - the size of the matrix.

A=zeros(n,n);

if nargin==2;
    m=1;
    for i=1:n
        A(i,i)=i;
        for j=2:n
            A(j,j-1)=m;
            A(j-1,j)=m;
        end
    end
end

if nargin==1
    for i=1:n
        A(i,i)=4;
        for j=2:n
            A(j,j-1)=1;
            A(j-1,j)=1;
        end
    end
end

end

