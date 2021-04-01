function T = fredkin(n)
%FREDKIN Create Fredkin matrix
%   F=FREDKIN(N) Create 3*N-by-3*N-bit Fredkin matrix.

if n==1
    T = blkdiag(eye(5),notm(1),1);
else
    T1 = fredkin(1);
    T=1;
    for i=1:n
        T=kron(T,T1);    
    end
end

