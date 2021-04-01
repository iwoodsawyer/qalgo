n = 1; % number of n input/output bits

% x is m-bits integer and output n-bits integer
% f = @(x,n) 0; % constant
%f = @(x,n) 1; % constant
f = @(x,n) x; % balanced (only for n=1)
%f = @(x,n) bitxor(x,2^n-1); % balanced

const = deutsch(f,n)
