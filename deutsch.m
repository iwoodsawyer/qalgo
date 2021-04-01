function const = deutsch(f,n) 
%DEUTSCH Deutsch algorithm
%  CONST = DEUTCH(F,N) determines if the function F is constant or
%  balanced, where F must be of the form F(X,N) where X is a N-bits
%  integer, and N is numbers of input and output bits.

if nargin == 1
    n = 1;
end

psi = kron(dec2vec(0,n),dec2vec(1,n));
Uf = ufm(f, n, n);
H = hadamard(n);
I = identity(n);
%psi = (kron(H,I)*(Uf*(kron(H,H)*psi)));
psi = Uf*kronmult({H,H},psi);
psi = kronmult({H,I},psi);
[psi,dec] = measure_subspace(psi, 1:n);
const = ~dec;



