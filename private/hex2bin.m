function bin=hex2bin(hex,n)
%HEC2BIN Convert hexidecimal string to a binary string.
%   HEX2BIN(D) returns the binary representation of H as a string.
%
%   HEX2BIN(D,N) produces a binary representation with at least N bits.
%
%   Example
%      hex2bin('1A') returns '11010'
%
%   See also HEX2DEC, DEC2BIN.

d = hex2dec(hex);
if ((nargin < 2) || isempty(n))
    bin = dec2bin(d);
else
    bin = dec2bin(d,n);
end