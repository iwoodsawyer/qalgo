function I = identity(n)
%IDENTITY Create identity matrix
%   I=IDENTITY(N) Create N-bit identity matrix.

I = speye(2^n);