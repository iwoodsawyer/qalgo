function phi = bin2vec(bin)
%BIN2VEC Converts a single (non-superposed) binary state description of a
%        register to a vector state representation.
%   PHI=BIN2VEC(B) interprets the binary string B and returns in PHI the
%   equivalent vector state representation. 
%
%   Example
%       phi=bin2vec('011')
%
%   See also BIN2DEC, DEC2VEC.

dec=bin2dec(bin);
phi=dec2vec(dec,length(bin));