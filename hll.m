function [x,p,a] = hll(A,b,n,r,max_meas) 
%HHL Harrow, Hassidim, and Lloyd linear system solver algorithm
%  [X,P,A] = HLL(A,B,N,R) determines the probability P of the solution 
%  X=A/B, where A must be N-by-N and B must be N-x-1.
%  V are the eigenvalues of matrix A.
%  R is scale factor for ancilla rotation.
%  A is the ancilla bit for the correctness of measurement.

%  Algorithm can currenly deal only with eigenvalues of value 2^k 

% input checks
if nargin < 5
    max_meas = 1e6;
end
assert(max_meas > 1);
if nargin < 4
    r = 0;
end
assert(r >= 0);

nl = log2(size(A,1));
nm = log2(size(A,2));
assert(nl==floor(nl));
assert(nl==nm, 'A is not square');
assert(abs(sum(norm(A*A'))-1)<sqrt(eps), 'A is not unitary');

if nargin < 3
    n = size(A,1);
end
assert(n>0,'N is not > 0');

if nargin < 2
    b = dec2vec(1,nl);
end
assert(iscolumn(b),'B is not column vector');
nb = log2(size(b,1));
assert(nb==floor(nb));
assert(nl==nb);
assert(abs(sum(norm(b))-1)<sqrt(eps), 'B is not normalised');


% initialize state
psi = kron(kron(b,dec2vec(0,n)),dec2vec(0,1));

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
    %psi = kron(Uctlm,I)*psi + psi;
    Uctlm = uctlm_mult(A, k, n, nb);
    psi = kronmult([Uctlm; {I}],psi) + psi;
end

% Inverse Quantum Fourier Transform (IQFT(
%Q = qft(n);
%psi = kron(kron(Ib,Q'),I)*psi;
%psi = kronmult([{Ib}; Q'; {I}],psi);
p = 2^nb;
l = 2^n;
m = l*2;
for i = 1:p
    for j = 1:2
        psi((i-1)*m+j:2:i*m) = fft(psi((i-1)*m+j:2:i*m))./sqrt(l);
    end
end

% eigenvalue inversion (simplified swap for now)
for k = 1:floor(n/2)
    %S = swap(k,n-k+1,n);
    %psi = kron(kron(Ib,S),I)*psi;
    [SI,SX,SY,SZ] = swap_mult(k,n-k+1,n);
    psi = kronmult([{Ib}; SI; {I}],psi) +...
          kronmult([{Ib}; SX; {I}],psi) +...
          kronmult([{Ib}; SY; {I}],psi) +...
          kronmult([{Ib}; SZ; {I}],psi);
end

% ancilla Rotation 
c = pi/((2^(r+2)));
for k = 1:n
    theta = c*2^(1-k);
    %Uinv = uinv(theta, k, n);
    %psi = kron(Ib,Uinv)*psi + psi;
    Uinv = uinv_mult(theta, k, n);
    psi = kronmult([{Ib}; Uinv],psi) + psi;
end

% undo eigenvalue inversion (simplified swap for now)
for k = floor(n/2):-1:1
    %S = swap(k,n-k+1,n);
    %psi = kron(kron(Ib,S),I)*psi;
    [SI,SX,SY,SZ] = swap_mult(k,n-k+1,n);
    psi = kronmult([{Ib}; SI; {I}],psi) +...
          kronmult([{Ib}; SX; {I}],psi) +...
          kronmult([{Ib}; SY; {I}],psi) +...
          kronmult([{Ib}; SZ; {I}],psi);
end

% Quantum Fourier Transform (QFT)
%psi = kron(kron(Ib,Q),I)*psi;
%psi = kronmult([{Ib}; Q; {I}],psi);
for i = 1:p
    for j = 1:2
        psi((i-1)*m+j:2:i*m) = ifft(psi((i-1)*m+j:2:i*m)).*sqrt(l);
    end
end

% undo control gates
for k = n:-1:1
    %Uctlm = uctlm(A', k, n, nb);
    %psi = kron(Uctlm,I)*psi + psi;
    Uctlm = uctlm_mult(A', k, n, nb);
    psi = kronmult([Uctlm; {I}],psi) + psi;
end

% undo Hadamard
%psi = kron(kron(Ib,H),I)*psi;
psi = kronmult([{Ib}; H; {I}],psi);

% measurement until a=1
a = 0;
k = 1;
while (abs(a) < sqrt(eps) && k<max_meas)
    [phi,a] = measure_subspace(psi,n+nb+1);
    k = k+1;
end
x = kron(kron(Ib,dec2vec(0,n)'),dec2vec(1,1)')*phi;
x = normalize(x);
p = abs(x).^2;
x = real(sign(x).*abs(x));
end


