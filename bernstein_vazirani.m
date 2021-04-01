function dec = bernstein_vazirani(f,m) 
%BERNSTEIN_VAZIRANI Bernstein-Vazirani algorithm
%  S = BERNSTEIN_VAZIRANI(F,M) finds the secret S in the function 
%  F = bindot(X,S), where S is not the zero-vector, where F must be of the
%  form F(X,N) with X is a M-bits integer, and M is numbers of input bits.

n = 1;
if nargin < 2
    m = 1;
end

psi = kron(dec2vec(0,m),dec2vec(1,n));
Uf  = ufm(f, m, n);
Hm  = hadamard(m);
Hn  = hadamard(n);
In  = identity(n);
psi = (kron(Hm,In)*(Uf*(kron(Hm,Hn)*psi)));
[psi,dec] = measure_subspace(psi, 1:m);











