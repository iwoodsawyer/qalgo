function [lambda,phi,a] = pe(A,psi,d,p)
%PE Phase Estimation algorithm
%[LAMBDA,PHI,A] = PE(A,PSI,D,P) determines the eigenvalue of A and initial
%  phase PSI, where A must be N-by-N and PSI must be N-x-1.
%  D is number of digits accuracy, and P is the probability of failure.
%  LAMBDA is the measured eigenvalue, and PHI is the measured top register.
%  A is the ancilla bit for the correctness of measurement.

% input checks
nl = log2(size(A,1));
nm = log2(size(A,2));
assert(nl==floor(nl));
assert(nl==nm, 'A is not square');
assert(abs(sum(norm(A*A'))-1)<sqrt(eps), 'A is not unitary');

if nargin < 4
    p = 1;
end
assert(p<=1,'P is not <= 1');

if nargin < 3
    d = size(A,1);
end
assert(d>0,'D is not > 0');

if nargin < 2
    psi = dec2vec(1,nl);
end
assert(iscolumn(psi),'PSI is not column vector');
nb = log2(size(psi,1));
assert(nb==floor(nb));
assert(nl==nb);
assert(abs(sum(norm(psi))-1)<sqrt(eps), 'PSI is not normalised');

% add qubits to reduce the probability of failure for top register
n = d + ceil(log2(1/2 + 1/(2*p)));

% initialize state
phi = kron(psi,dec2vec(0,n));

% apply Hadamard
Ib = identity(nb);
% H = hadamard(n);
% phi = kron(Ib,H)*phi;
H = hadamard_mult(n);
phi = kronmult([{Ib}; H],phi);

% apply control gates
for k = 1:n
    %Uctlm = uctlm(A, k, n, nb);
    %phi = Uctlm*phi + phi;
    Uctlm = uctlm_mult(A, k, n, nb);
    phi = kronmult(Uctlm,phi) + phi;
end

% apply inverse Quantum Fourier Transform
% Q = qft(n);
% phi = kron(Ib,Q')*phi;
% phi = kronmult([{Ib}; Q'],phi);
p = 2^nb;
l = 2^n;
for j = 1:p
    phi((j-1)*l+1:j*l) = fft(phi((j-1)*l+1:j*l))./sqrt(l);
end



% measure
[~,a] = measure_subspace(phi,1:nb);
[phi,lambda] = measure_subspace(phi,nb+1:nb+n);

% convert
lambda = exp(2*pi*1i*lambda/(2^n));

end