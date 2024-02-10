%% 2x2 PauliZ
A = pauliz(1);
psi = dec2vec(1,1);

p = 1; % probability of failure
n = 2;
m = 4; % number of digits accuracy

d = eig(A)
d1 = pe(A,psi,n,p)
d2 = pe_itr(A,psi,m)


%% 2x2 Phase
A = phase(1);
psi = dec2vec(1,1);

d = eig(A)
d1 = pe(A,psi,n,p)
d2 = pe_itr(A,psi,m)


%% 2x2 Matrix
A = [3, -1; -1, 3]; 
%psi = pauliz(1)*hadamard(1)*dec2vec(0,1)
psi = hadamard(1)*dec2vec(0,1)
[U,d,n,t] = uexpm(A);

d

v1 = pe(U,psi,n,p)
d1 = (log(abs(v1))+ 1i*(angle(v1) + 2*pi*(angle(v1) < 0)))./(t*1i)

v2 = pe_itr(U,psi,m)
d2 = (log(abs(v2))+ 1i*(angle(v2) + 2*pi*(angle(v2) < 0)))./(t*1i)


%% 4x4 Matrix
A = [15 9 5 -3; 9 15 3 -5; 5 3 15 -9; -3 -5 -9 15]./4;
psi = [1; 1; -1; 1]./2;
[U,d,n,t] = uexpm(A);

d

v1 = pe(U,psi,n,p)
d1 = (log(abs(v1))+ 1i*(angle(v1) + 2*pi*(angle(v1) < 0)))./(t*1i)

v2 = pe_itr(U,psi,m)
d2 = (log(abs(v2))+ 1i*(angle(v2) + 2*pi*(angle(v2) < 0)))./(t*1i)
