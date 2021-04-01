function phi=dec2vec(dec,n)
%DEC2VEC Convert single decimal state to vector representation
%   PHI=DEC2VEC(D,N) interprets the single decimal D and returns in PHI the
%   equivalent vector state representation. The N is max bits.
%
%   Example
%       phi=dec2vec(3,4)
%
%   See also BIN2DEC, BIN2VEC.

assert(dec<=2^n-1, 'Decimal too big for register size n!');

phi=zeros(2^n,1);
phi(dec+1)=1; %compensate for matlab indexing with +1