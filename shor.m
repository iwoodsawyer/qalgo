function [p,q] = shor(N) 
%SHOR Shor's factoring algorithm
%  [P,Q]=SHOR(N) factorizes a positive integer N into N=P*Q if it exists. 

assert(N==floor(N),'N is not an integer');
assert(N < 255,'N=255 is the biggest number');
assert(~ispn(N),'N is prime number');
assert(~ispnp(N),'N is prime power');

n = ceil(log2(N));
m = 2*n;

k = 1;
found = false;
while (~found && k<(N-2))
    % randomly choose integer with 1<a<N
    a = floor(2+(N-3)*rand);
    
    % check if a is not 1
    p = gcd(int64(a),int64(N));
    %p = gcd(a,N);
    if p~=1
       found = true; 
       q = N/p; 
       continue;
    end
    
    % modulo operation
    f = @(x,y) mod((y*a^x),N);
    
    % initilize state vector
    psi = kron(dec2vec(0,m),dec2vec(0,n));
    
    % apply hadamard
    %Hm  = hadamard_mult(m);
    Hm  = hadamard(m);
    In  = identity(n);
    %psi = kronmult({Hm{:} In},psi);
    psi = kronmult({Hm,In},psi); 
    
    % apply operator
    %Uf  = ufm(f, m, n);
    Uf  = ufam(f, m, n);
    psi = Uf*psi; 
    
    % measure
    psi = measure_subspace(psi, m+1:m+n);
    
    % apply QFT
    Qm  = qft(m);
    psi = kronmult({Qm,In},psi);
    %psi = kron(Qm,In)*psi;
    
    % measure
    [psi,x] = measure_subspace(psi, 1:m);
    
    % determine period
    if (x==0)
        continue; % try agian, no information about the period
    else
        r = int64(cfa(x,2^m));
    end
    
    % determine factors
    if mod(r,2) % is odd
       continue; % try agian, period must be even
    else % is even
        d = int64(a)^(r/2);
        if d==mod(-1,N)
            continue; % try agian, no information about N
        else
            p = gcd(d+1,int64(N));
            if ((p > 1) && (p < N))
                found = true;
                q = N/p;
            else
                p = gcd(d-1,int64(N));
                if ((p > 1) && (p < N))
                    found = true;
                    q = N/p;
                else
                    continue; % try agian, no factor other then N and 1 obtained
                end 
            end
        end
    end
end

if ~found
    error('Period not found.')
end
end

function err=ispn(n)
%ISPN Is N a prime number
err=true;
if n<=1
    return
end
for i=2:floor(sqrt(n))
    if mod(n,i)==0
        err=false;
        return
    end
end
end

function err=ispnp(n)
%IPNP Is N is a prime power
err=false;
i=2;
f=0;
while (i<=floor(sqrt(n)) && f==0)
    if mod(n,i)==0
        f=i;
    end
    i=i+1;
end
i=2;
while (f^i<=n && ~err)
    if f^i==n
        err = true;
    end
    i=i+1;
end
end

function r = cfa(x,qmax)
%CFA Continued Fraction Algorithm
r=1;
i=0;
while (r>0)
    i=i+1;
    lambda(i)=floor(qmax/x);
    r=qmax-lambda(i)*x;
    qmax=x;
    x=r;
end

for i=length(lambda):-1:1
    a=lambda(i);
    b=1;
    for j=i-1:-1:1
        c=a;
        a=lambda(j)*a+b;
        b=c;
    end
    r=max(r,a);
end
end




