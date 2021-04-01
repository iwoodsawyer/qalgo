function H = hadamard_mult(n)
%HADAMARD Create Hadamard matrix
%   H=HADAMARD(N) Create N-bit Hadamard matrix.

if n==1
    H = [ 1 1; 1 -1]/sqrt(2);
else
    H=cell(n,1);
    for i=1:n
        H{i} = hadamard(1);    
    end
end

