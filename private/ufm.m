function Uf = ufm(f,m,n)
%UFM Make unitary UF matrix from function F
%   UF = UFM(F,M,N) makes unitary UF from function F, where F must be of
%   the form F(X,N) where X is a M-bits integer, and M and N is the numbers
%   of input and output bits.

% i is integer version of string x
% j is integer version of string b
% f:i -> R^n

k = 2^m;
l = 2^n;
s = k*l;
ins=zeros(s,1);
outs=zeros(s,1);

for i=0:k-1
    fi = feval(f,i,n);
    ins(i*l+(1:l))  = bitshift(i,n) + (0:l-1) + 1;
    outs(i*l+(1:l)) = bitshift(i,n) + bitxor((0:l-1), fi) + 1;
end

Uf = full(sparse(ins,outs,ones(s,1),s,s));


