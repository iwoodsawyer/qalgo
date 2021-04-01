function C = binmult(A,B)
%BINNUMT Multiplication of binary matrices
%  C = BINMULT(A,B) multiplies binary matrcied A*B=C.

[ma,na] = size(A);
[mb,nb] = size(B);
assert(na==mb)

C = char(zeros(ma,nb));
for i = 1:ma
    for j = 1:nb
        for k = 1:na
            C(i,j) = dec2bin(bitxor(bin2dec(C(i,j)),(bin2dec(A(i,k))*bin2dec(B(k,j)))));
        end
    end
end

