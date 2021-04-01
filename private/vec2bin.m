function bin = vec2bin(phi)
%VEC2BIN Convert vector representation to single binary state description
%   B=VEC2BIN(PHI) interprets the equivalent vector state representation
%   PHI and returns single (non-superposed) binary state description B.

n=log2(size(phi,1));
assert(n==floor(n));

in=find(phi~=0); %nonzero
dec=0;
for i=1:length(in)
    dec = dec+in(i)-1;
end
bin = dec2bin(dec,n);
