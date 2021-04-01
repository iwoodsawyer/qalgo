function phi = num2vec(num)
%HEX2VEC Converts a decimal number to a vector state representation.
%   PHI=NUM2VEC(F) interprets the decimal number F and returns in PHI
%   the equivalent vector state representation.
%
%   Example
%       phi=num2vec(single(pi))
%
%   See also BIN2VEC, DEC2VEC.

hex=num2hex(num);
bin=hex2bin(hex);
dec=bin2dec(bin);
phi=dec2vec(dec,length(bin));