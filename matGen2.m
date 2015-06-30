function[A]=matGen2(n,m,k)
%This function generates the matrix type as described in Example 2

A=zeros(n,n);

for i=1:n
    A(i,i)=m+i;
    for j=2:n
        r=mod(j,2);
        if r==0
            A(j,j-1)=sqrt(j-(.5*j));
            A(j-1,j)=sqrt(j-(.5*j));
        else
            A(j,j-1)=k*sqrt(m+(j-1)-(.5*(j-1)));
            A(j-1,j)=k*sqrt(m+(j-1)-(.5*(j-1)));
        end
    end
end

end

