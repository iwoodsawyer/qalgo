n = 3; % number of n input/output bits
c = 5;
f = @(x) (x==c); % balanced
c = grover(f,n) % c=5


n = 3; % number of n input/output bits
c1 = 5;
c2 = 3;
f = @(x) (x==c1 || x==c2); % balanced
c = grover(f,n,2)


n = 4; % number of n input/output bits
c1 = 5;
c2 = 3;
c3 = 13;
f = @(x) (x==c1 || x==c2 || x==c3); % balanced
c = grover_gen(f,n)


n = 4; % number of n input/output bits
min = 5;
f = @(x) ceil(((x-min)^2)/8);
c = grover_min(f,n)
