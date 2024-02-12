%% 2x2 Matrix with 2^k eigenvalues
A = [3, -1; -1, 3]; 
%b = pauliz(1)*hadamard(1)*dec2vec(0,1)
%b = hadamard(1)*dec2vec(0,1)
b = dec2vec(0,1);
[U,d,n,t] = uexpm(A);

x = A\b;
x = normalize(x)
p = abs(x).^2 % probability 

r = 5;
[x1,p1] = hll(U,b,n,r)

%% 2x2 Matrix with (2^k)/3 eigenvalues
A = [1, -1/3; -1/3, 1];
%b = pauliz(1)*hadamard(1)*dec2vec(0,1)
%b = hadamard(1)*dec2vec(0,1)
b = dec2vec(0,1);
[U,d,n,t] = uexpm(A);

x = A\b;
x = normalize(x)
p = abs(x).^2 % probability 

r = 5;
[x1,p1] = hll(U,b,n,r)

%% 2x2 Matrix with (2^k)+1 eigenvalues
A = [4, 1; 1, 4];
%b = pauliz(1)*hadamard(1)*dec2vec(0,1)
%b = hadamard(1)*dec2vec(0,1)
b = dec2vec(0,1);
[U,d,n,t] = uexpm(A);

x = A\b;
x = normalize(x)
p = abs(x).^2 % probability 

r = 5;
[x1,p1] = hll(U,b,n,r)

%% 4x4 Matrix with 2^k eigenvalues
A = [15 9 5 -3; 9 15 3 -5; 5 3 15 -9; -3 -5 -9 15]./4;
b = [1; 1; 1; 1]./2;
[U,d,n,t] = uexpm(A);

x = A\b;
x = normalize(x)
p = abs(x).^2 % probability 

r = 5;
[x1,p1] = hll(U,b,n,r)
