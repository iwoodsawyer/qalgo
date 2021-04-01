m = 4; % number of m input bits
n = 1; % number of n output bits

% x is m-bits integer and output n-bits integer
%f = @(x,n) 1; % constant
%f = @(x,n) 0; % constant
%f = @(x,n) ~mod(x,2); % balanced
f = @(x,n) ~sign(x); % (if n=1 balanced, for n>1 nor balanced nor constant, for high n high probality constant)
Uf  = full(ufm(f, m, n))

[const,bal] = deutsch_jozsa(f,m,n)
