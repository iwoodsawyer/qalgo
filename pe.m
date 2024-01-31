function [lambda,phi] = pe(A,psi,p)
%PE Phase Estimation algorithm
%[LAMBDA,PHI] = PE(A,PSI,P) determines the eigenvalue of A and initial
%  phase PSI, where A must be N by N-bits and PSI must be 1 x N-bits.
%  P is the probability of failure of this algorithm.

[m,n] = size(A);
if m~=n
   error('A must be square'); 
end
if nargin < 3
    p = 1;
end
if nargin < 2
    psi = zeros(n,1);
    psi(end) = 1;
end

% calculate the number of qubits for top register
n = n + ceil(log2(2 + 1/(2*p)));

% apply hadamard
phi = dec2vec(0,n);
% H = hadamard(n);
% phi = H*phi;
H = hadamard_mult(n);
phi = kronmult(H,phi);

% apply control U^(2^(k-1))
for k = 1:n
    %phi = uctl(A, psi, k, n) * phi;
    phi = kronmult(uctl(A, psi, k, n),phi);
end

% apply inversed Quantum Fourier Transform
%phi = qft(n)'*phi;
phi = fft(phi)./sqrt(length(phi));

[phi,lambda] = measure(phi);
lambda = exp(2*pi*1i*lambda/(2^n));

end