function [x,s] = hll(A,b,v) 
%HHL Harrow, Hassidim, and Lloyd linear system solver algorithm
%  X = HLL(A,B,V)

[m,n] = size(A);
nb = log2(length(b));
if m~=n
   error('A must be square'); 
end

% State Preperation
psi = kron(b,kron(dec2vec(0,n),dec2vec(0,1)));


%% IQFT

% apply Hadamard
I = identity(1);
Ib = identity(nb);
H = hadamard(n);
psi = kron(kron(Ib,H),I)*psi;

% apply control gates
for k = 1:n
    Uctlm = uctlm2(A, k, n, nb);
    psi = kron(Uctlm,I)*psi;
end

% apply inverse Quantum Fourier Transform
psi = kron(kron(Ib,qft(n)'),I)*psi;


%% Ancilla Rotation 
for k = 1:length(v)
    Uinv = uinv(v, k, n);
    psi = kron(Ib,Uinv)*psi;
end

%% QFT

% Quantum Fourier Transform
psi = kron(kron(Ib,qft(n)),I)*psi;

% apply control gates
for k = n:-1:1
    Uctlm = uctlm2(A', k, n, nb);
    psi = kron(Uctlm,I)*psi;
end

% apply Hadamard
psi = kron(kron(Ib,H),I)*psi;

% Measurement
[phi,a] = measure_subspace(psi,n+nb+1);

if a > sqrt(eps)
    disp('WELL');
    s = vec2struct(phi);
    x = [s.alpha].*conj([s.alpha]);
else
    disp('NOTHING');
    s = 0;
    x = 0;
end

end


