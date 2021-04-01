function T = toffoli(n)
%TOFFOLI Create Toffoli matrix
%   T=TOFFOLI(N) Create 3*N-by-3*N-bit Toffoli matrix.

if n==1
    T = blkdiag(eye(6),notmat(1));
else
    T1 = toffoli(1);
    T=1;
    for i=1:n
        T=kron(T,T1);    
    end
end

