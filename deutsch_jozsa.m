function [const,bal] = deutsch_jozsa(f,m,n) 
%DEUTSCH_JOZSA Deutsch-Jozsa algorithm
%  [CONST,BAL] = DEUTCH_JOZSA(F,M,N) determines if the function F is
%  constant or balanced, where F must be of the form F(X,N) where X is a
%  M-bits integer, and M and N is numbers of input and output bits.

if nargin < 3
    n = 1;
end
if nargin < 2
    m = 1;
end

psi = kron(dec2vec(0,m),dec2vec(1,n));
Uf  = ufm(f, m, n);
Hm  = hadamard(m);
Hn  = hadamard(n);
In  = identity(n);
%psi = (kron(Hm,In)*(Uf*(kron(Hm,Hn)*psi)));
psi = Uf*kronmult({Hm,Hn},psi);
psi = kronmult({Hm,In},psi);
[psi,dec] = measure_subspace(psi, 1:m);
const = ~dec;
bal = (dec==1);












