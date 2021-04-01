function Uf = ufam(f,m,n)
%UFAM Make unitary UF matrix from function F = A^X MOD N
%   UF = UFAM(F,M,Y) makes unitary UF from function F = A^X MOD N, where F
%   must be of the form F(X,Y) where X is a M-bits integer. Stable for
%   large values of X, because of F(i) = A^X(i) MOD N = F(i-1)A MOD N.

k = 2^m;
l = 2^n;
s = k*l;
ins=zeros(s,1);
outs=zeros(s,1);

fi = zeros(k,1);
for i=0:k-1
    % trick y(i) = a^x(i) mod N = y(i-1)*a mod N
    if i==0
        fi(i+1) = feval(f,i,1);
    else
        fi(i+1) = feval(f,1,fi(i));
    end
    ins(i*l+(1:l))  = bitshift(i,n) + (0:l-1) + 1;
    outs(i*l+(1:l)) = bitshift(i,n) + bitxor((0:l-1), fi(i+1)) + 1;
end

Uf = sparse(ins,outs,ones(s,1),s,s);


