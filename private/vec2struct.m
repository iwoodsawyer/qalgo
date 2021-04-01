function s = vec2struct(phi)

n=log2(size(phi,1));
assert(n==floor(n));

in=find(abs(phi) > sqrt(eps)); %nonzero

l = length(in);
s(l).bin='';
s(l).alpha=1;
for i=1:l
    s(i).bin = dec2bin(in(i)-1, n);
    s(i).alpha = phi(in(i));
end

