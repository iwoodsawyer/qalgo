function [x,a] = hll(A,b,n,theta) 
%HHL Harrow, Hassidim, and Lloyd linear system solver algorithm
%  [X,A] = HLL(A,B,V) determines the probability of the solution X=A/B 
%  where A must be N by N-bits and full rank and B must be 1 x N-bits.
%  V are the eigenvalues of matrix A.

% input checks
[m,nl] = size(A);
assert(m==nl, 'A not aquare.');
assert(abs(sum(norm(A))-1)<sqrt(eps), 'A not normalised.');
b = normalize(b);
nb = log2(size(b,1));
assert(nb==floor(nb));
assert(abs(sum(norm(b))-1)<sqrt(eps), 'B not normalised.');


% add qubits to reduce the probability of failure for top register
p = 1;
n = n + ceil(log2(1/2 + 1/(2*p)));

% initialize state
psi = kron(kron(b,dec2vec(0,n)),dec2vec(0,1));


%% IQFT

% apply Hadamard
I = identity(1);
Ib = identity(nb);
H = hadamard(n);
psi = kron(kron(Ib,H),I)*psi;
% H = hadamard_mult(n);
% psi = kronmult([{Ib}; H; {I}],psi);

% apply control gates
for k = 1:n
    Uctlm = uctlm(A, k, n, nb);
    psi = kron(Uctlm,I)*psi + psi;
    %Uctlm = uctlm_mult(A, k, n, nb);
    %psi = kronmult([Uctlm; {I}],psi) + psi;
end

% apply inverse Quantum Fourier Transform
Q = qft(n);
psi = kron(kron(Ib,Q'),I)*psi;
%psi = kronmult([{Ib}; Q'; {I}],psi);
% p = 2^nb;
% l = 2^n;
% m = p*l;
% for i = 1:p
%     for j = 1:2
%         psi((i-1)*m+j:2:i*m) = fft(psi((i-1)*m+j:2:i*m))./sqrt(l);
%     end
% end

%% Eigenvalue inversion (simplified swap for now)
% for k = 1:floor(n/2)
%     S = swap(k,n-k+1,n);
%     psi = kron(kron(Ib,S),I)*psi;
% end

%% Ancilla Rotation 
for k = 1:length(theta)
    Uinv = uinv(theta(k), k, n);
    psi = kron(Ib,Uinv)*psi + psi;
    %Uinv = uinv_mult(theta(k), k, n);
    %psi = kronmult([{Ib}; Uinv],psi) + psi;
end

%% QFT

% Quantum Fourier Transform
psi = kron(kron(Ib,Q),I)*psi;
%psi = kronmult([{Ib}; Q; {I}],psi);
% for i = 1:p
%     for j = 1:2
%         psi((i-1)*m+j:2:i*m) = ifft(psi((i-1)*m+j:2:i*m)).*sqrt(l);
%     end
% end

% apply control gates
for k = n:-1:1
    Uctlm = uctlm(A', k, n, nb);
    psi = kron(Uctlm,I)*psi + psi;
    %Uctlm = uctlm_mult(A', k, n, nb);
    %psi = kronmult([Uctlm; {I}],psi) + psi;
end

% apply Hadamard
psi = kron(kron(Ib,H),I)*psi;
%psi = kronmult([{Ib}; H; {I}],psi);

% Measurement
[~,a] = measure_subspace(psi,n+nb+1);
[~,~,~,x] = measure_subspace(psi,1:nb);
x = normalize(x);
x = x.^2;

end


