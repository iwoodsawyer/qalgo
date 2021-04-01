function phi = hex2vec(hex)
%HEX2VEC Converts a hexidecimal string to a vector state representation.
%   PHI=HEX2VEC(H) interprets the hexidecimal string H and returns in PHI
%   the equivalent vector state representation.
%
%   Example
%       phi=hex2vec('1A')
%
%   See also BIN2VEC, DEC2VEC.

bin=hex2bin(hex);
dec=bin2dec(bin);
phi=dec2vec(dec,length(bin));