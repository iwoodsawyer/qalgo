function c = grover(f,n,p) 
%GROVER Grover's search algorithm
%  C = GROVER(F,N) for searching an unsorted database F for entry C, where
%  F must be of the form F(X) with X is a N-bits integer. The default for N
%  is 1. The function F returns 1 for X=C and 0 for other X. The intger P
%  is the number of occurances of C in function F. The default for P is 1.

if nargin <3
    p = 1;
end
if nargin <2
    n = 1;
end
maxiter = floor((pi/4)*sqrt(2^n / p));

phi = dec2vec(0,n);
H = hadamard(n);
D = iam(n);
Ui = uim(f,n);
phi = H*phi;
for iter=1:maxiter
    phi=D*(Ui*phi);
end
phi = measure(phi);
c = vec2dec(phi);
