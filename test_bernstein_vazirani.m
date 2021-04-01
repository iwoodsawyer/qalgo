m = 3; % number of m input bits
s = 7;

% x is m-bits integer and output 1-bits integer
f = @(x,n) bin2dec(bindot(dec2bin(x,m),dec2bin(s,m))); 


Uf  = full( ufm(f, m, 1) )

s =  bernstein_vazirani(f,m)
