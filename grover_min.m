function c = grover_min(f,n) 
%GROVER Generic Grover's minimum search algorithm
%  C = GROVER(F,N) for searching minimum unsorted database F for entry C,
%  where F must be of the form F(X) with X is a N-bits integer. The default
%  for N is 1. The function F returns 1 for X=C and 0 for other X.

if nargin <2
    n = 1;
end
N = 2^n;
H = hadamard(n);
D = iam(n);

y = N-1;
c = 0;
m = 1;
lambda = 6/5;
found = false;
while ~found
    j = ceil(floor(m)*rand);
    g = @(x,y) f(x)<y;
    
    phi = dec2vec(0,n);
    phi = H*phi;
    Ui = uim(g,n,y);
    for iter=1:j
        phi=D*(Ui*phi);
    end
    [phi,c1] = measure(phi);
    
    y1 = f(c1); 
    if y1 < y
        y = y1;
        c = c1;
        m = 1;
    elseif m >= sqrt(N);
        found = true;
    else
        m = min(lambda*m, sqrt(N));
    end
end
