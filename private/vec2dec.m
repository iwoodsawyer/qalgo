function dec = vec2dec(phi)
%VEC2DEC Convert vector representation to single decimal state
%   D=VEC2DEC(PHI) interprets the equivalent vector state representation
%   PHI and returns the single decimal D.

n=log2(size(phi,1));
assert(n==floor(n));

in=find(phi~=0); %nonzero
dec=0;
for i=1:length(in)
    dec = dec+in(i)-1;
end

