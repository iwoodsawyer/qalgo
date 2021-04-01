function H = xorm(n)
%XORM Create XOR matrix
%   H=XORM(N) Create N-by-2N-bit XOR matrix.

if n==1
    H = [ 1 0 0 1; 0 1 1 0];
else
    H1 = xorm(1);
    H=1;
    for i=1:n
        H=kron(H,H1);    
    end
end

