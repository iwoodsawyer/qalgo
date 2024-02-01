function [x,s] = hll(A,b,v) 
%HHL Harrow, Hassidim, and Lloyd linear system solver algorithm
%  X = HLL(A,B,V) determines the probability of the solution X=A/B where A 
%  must be N by N-bits and full rank and B must be 1 x N-bits.
%  V are the eigenvalues of matrix A.

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
%H = hadamard(n);
%psi = kron(kron(Ib,H),I)*psi;
H = hadamard_mult(n);
psi = kronmult([{Ib}; H; {I}],psi);

% apply control gates
for k = 1:n
      %Uctlm = uctlm(A, k, n, nb);
      %psi = kron(Uctlm,I)*psi;
      [Im,Upow] = uctlm_mult(A, k, n, nb);
      psi = kronmult([Im; {I}],psi) + kronmult([Upow; {I}],psi);
end

% apply inverse Quantum Fourier Transform
%Q = qft(n);
%psi = kron(kron(Ib,Q'),I)*psi;
%psi = kronmult([{Ib}; Q'; {I}],psi);
p = 2^nb;
l = 2^n;
m = p*l;
for i = 1:2
    for j = 1:p
        psi((i-1)*m+j:p:i*m) = fft(psi((i-1)*m+j:p:i*m))./sqrt(l);
    end
end

%% Ancilla Rotation 
v = v./min(v);
for k = 1:length(v)
    %Uinv = uinv(v, k, n);
    %psi = kron(Ib,Uinv)*psi;
    [Im,Uinv] = uinv_mult(v, k, n);
    psi = kronmult([{Ib}; Im],psi) + kronmult([{Ib}; Uinv],psi);
end

%% QFT

% Quantum Fourier Transform
%psi = kron(kron(Ib,Q),I)*psi;
%psi = kronmult([{Ib}; Q; {I}],psi);
for i = 1:2
    for j = 1:p
        psi((i-1)*m+j:p:i*m) = ifft(psi((i-1)*m+j:p:i*m)).*sqrt(l);
    end
end

% apply control gates
for k = n:-1:1
    %Uctlm = uctlm(A', k, n, nb);
    %psi = kron(Uctlm,I)*psi;
    [Im,Upow] = uctlm_mult(A', k, n, nb);
    psi = kronmult([Im; {I}],psi) + kronmult([Upow; {I}],psi);
end

% apply Hadamard
%psi = kron(kron(Ib,H),I)*psi;
psi = kronmult([{Ib}; H; {I}],psi);

% Measurement
[phi,a] = measure_subspace(psi,n+nb+1);

if a > sqrt(eps)
    disp('WELL');
    s = vec2struct(phi);
    x = ([s.alpha].*conj([s.alpha]))';
else
    disp('NOTHING');
    s = 0;
    x = 0;
end

end


