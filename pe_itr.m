function [lambda,phi] = pe_itr(A,psi,d)
%PE_ITR Iterative Phase Estimation algorithm
%[LAMBDA,PHI] = PE_ITR(A,PSI,D) determines the eigenvalue of A and initial
%  phase PSI, where A must be N-by-N and PSI must be 1-x-N.
%  D is number of digits accuracy.
%  LAMBDA is the measured eigenvalue, and PHI is the measured top register.

% input checks
assert(abs(sum(norm(A*A'))-1)<sqrt(eps), 'A is not unitary');
nb = log2(size(psi,1));
assert(nb==floor(nb));
if nargin < 3
    d = size(A,1);
end
assert(d>0,'D is not > 0');
if nargin < 2
    psi = dec2vec(1,nb);
end
assert(abs(sum(norm(psi))-1)<sqrt(eps), 'PSI is not normalised');

n = 1;
I = identity(1);
Ib = identity(nb);
H = hadamard(n);

lambda = 0;
num_meas = 1000;
for k = d:-1:1   
    % inital state
    phi = kron(kron(psi,dec2vec(0,n)),dec2vec(0,1)); % c q a
    
    % apply hadamard
    phi = kron(kron(Ib,H),I)*phi;
    
    % apply control gate
    Uctl = uctl1(A, k, nb);
    phi = kron(Uctl,I)*phi + phi;
    
    % apply rotation
    Urot = urot(-lambda*(2^(k-1)), 1, n);
    phi = kron(Ib,Urot)*phi + phi;
    
    % apply hadamard
    phi = kron(kron(Ib,H),I)*phi;
    
    % measurements
    Theta = zeros(num_meas,1);
    for j = 1:num_meas     
        [~,theta] = measure_subspace(phi,nb+1:nb+n);
        Theta(j) = round(theta);
    end
    
    % select frequent measurement
    theta = mode(Theta);
    lambda = lambda + theta/(2^k);
end

% convert
lambda = exp(2*pi*lambda*1i);
end