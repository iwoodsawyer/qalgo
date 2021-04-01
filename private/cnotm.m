function H = cnotm(n)
%CNOTM Create controlled NOT matrix
%   H=CNOTM(N) Create 2N-by-2N-bit CNOT matrix.

if n==1
    H = [ 1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0 ];
else
    H1 = cnotm(1);
    H=1;
    for i=1:n
        H=kron(H,H1);    
    end
end

