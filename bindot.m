function c = bindot(a,b)
%BINDOT Dot product in single binary state description
%   C=BINDOT(A,B) calculates dot product from single binary state
%   description A and B returns single binary state description C.

assert(length(a)==length(b));
c = 0;
for i=1:length(a)
    c = bitxor(c,bin2dec(a(i))*bin2dec(b(i)));
end
c = dec2bin(c,1);