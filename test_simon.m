n = 3; % number of n input/output bits
c = 5;

% x is m-bits integer and output n-bits integer
f = @(x,n) bitxor(x,c); % balanced

c = simon('func_simon',n) % c=5
c = simon2('func_simon',n) % c=5

c = simon('func_simon1',n) % c=6
c = simon2('func_simon1',n) % c=6

c = simon('func_simon2',n) % c=3
c = simon2('func_simon2',n) % c=3

c = simon('func_simon3',n) % c=0 (constant)
c = simon2('func_simon3',n) % c=0 (constant)

c = simon(f,n) % no period
c = simon2(f,n) % no period