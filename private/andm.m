function H = andm(n)
%ANDM Create AND matrix
%   H=ANDM(N) Create N-by-2N-bit AND matrix.

if n==1
    H = [ 1 1 1 0; 0 0 0 1];
else
    H1 = andm(1);
    H=1;
    for i=1:n
        H=kron(H,H1);    
    end
end

