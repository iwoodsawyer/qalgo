p = 0.25; % probability
m = 4; % number of dgits

A = pauliz(1);
psi = dec2vec(1,1);

d = eig(A)
d1 = pe(A,psi,p)
d2 = pe_itr(A,psi,m)


A = phase(1);
psi = dec2vec(1,1);

d = eig(A)
d1 = pe(A,psi,p)
d2 = pe_itr(A,psi,m)