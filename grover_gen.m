function c = grover_gen(f,n) 
%GROVER Generic Grover's search algorithm
%  C = GROVER(F,N) for searching an unsorted database F for entry C, where
%  F must be of the form F(X) with X is a N-bits integer. The default for N
%  is 1. The function F returns 1 for X=C and 0 for other X. 

if nargin <2
    n = 1;
end
N = sqrt(2^n);
H = hadamard(n);
D = iam(n);
Ui = uim(f,n);

m = 1;
lambda = 6/5;
found = false;
while ~found
    j = ceil(floor(m)*rand);
    
    phi = dec2vec(0,n);
    phi = H*phi;
    for iter=1:j
        phi=D*(Ui*phi);
    end
    [phi,c] = measure(phi);
    
    if f(c)
        found = true;
    else
        m = min(lambda*m, N);
    end
end
