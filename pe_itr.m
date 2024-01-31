function [lambda,phi] = pe_itr(A,psi,d)
%PE_ITR Iterative Phase Estimation algorithm
%[LAMBDA,PHI] = PE_ITR(A,PSI,D) determines the eigenvalue of A and initial
%  phase PSI, where A must be N by N-bits and PSI must be 1 x N-bits.
%  D is the number of digits accuracy of this algorithm.

[m,n] = size(A);
if m~=n
   error('A must be square'); 
end
if nargin < 3
    d = n;
end
if nargin < 2
    psi = zeros(n,1);
    psi(end) = 1;
end

% calculate the number of qubits for top register
n = 1;

lambda = 0;
for k = d:-1:1
    % apply hadamard
    phi = dec2vec(0,n);
    
    H = hadamard(n);
    phi = H*phi;
     
    % apply control U^(2^(k-1))
    Uctl = uctl1(A, psi, k);
    phi = Uctl*phi;
     
    % apply rotation
    Z = rotz(-lambda*(2^(k-1)),n);
    phi = Z*phi;
    
    % apply inversed Quantum Fourier Transform
    phi = H*phi;
    
    [phi,theta] = measure(phi);
    
    lambda = lambda + theta/(2^k);
end

lambda = exp(2*pi*1i*lambda);

end