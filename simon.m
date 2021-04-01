function c = simon(f,n) 
%SIMON Simon's periodicity algorithm
%  C = SIMON(F,N) find the periodic component C in the function F =
%  BITXOR(X,C), where F must be of the form F(X,N) with X is a N-bits
%  integer, and N is numbers of input and output bits.

if nargin <2
    n = 1;
end

s = char(zeros(0,n));
for i=1:4*n % repeating Simon’s algorithm no more than 4n times
    psi = kron(dec2vec(0,n),dec2vec(0,n));
    Hn  = hadamard(n);
    In  = identity(n);
    Uf  = ufm(f, n, n);
    %psi = (kron(Hn,In)*(Uf*(kron(Hn,In)*psi)));  
    psi = Uf*kronmult({Hn,In},psi);
    psi = kronmult({Hn,In},psi);
    %[~,~,bin] = measure_subspace(psi, 1:n);
    [psi,~,bin] = measure(psi);
    bin = bin(1:n);
    
    s = [s; bin]; % add solution
    [r,s]=binrank(s); % return rank and the row-reduced gauss-jordan form
end

if r==n % found n distinct rows
    warning('No periodic function found!')
    c = 0; % c = bin2dec(binsolve(s,dec2bin(zeros(n,1)))); (only 1 solution c=0)
elseif r==n-1 % found n-1 distinct rows
    % column pivoting on the matrix to gauss-jordan form
    idx = 1:n;
    for i = 1:r
        if ~bin2dec(s(i,i))
            for j = i+1:1:n+1;
                if bin2dec(s(i,j))
                    idx([i j]) = idx([j i]); 
                    t = s(:,i);
                    s(:,i) = s(:,j);
                    s(:,j) = t;
                    break
                end
            end
        end
    end
    
    % calculate the solution of the Homogeneous equations
    c = [s(1:r,r+1)' '1'];
    
    % reverse column pivoting
    c = bin2dec(c(idx));  
elseif r==0
    warning('Constant function found!')
    c = 0;
else
    error('Undetermined function. Maximum iterations reached! (4*n)')
end






