function H = hadamard(n)
%HADAMARD Create Hadamard matrix
%   H=HADAMARD(N) Create N-bit Hadamard matrix.

if n==1
    H = [ 1 1; 1 -1]/sqrt(2);
else
    H1 = hadamard(1);
    H=1;
    for i=1:n
        H=kron(H,H1);    
    end
end

