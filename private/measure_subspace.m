function [phi,dec,bin,pi] = measure_subspace(psi, bits)
%MEASURE_SUBSPACE Measurement on a subspace of qubits of register psi.
%   PHI=MEASURE_SUBSPACE(PSI,BITS) specify which bits to measure with a N
%   indices, where N is number of qubits in to be measured. For example,
%   MEASURE_SUBSPACE(PSI,[2 4]) will measure qubits 2 and 4 of a 4-qubit
%   register state in PSI. PHI is the collapsed final state after the
%   observation.
%
%   [PHI,DEC]=MEASURE_SUBSPACE(PSI,BITS) gives result in decimal state
%   description.
%
%   [PHI,DEC,BIN]=MEASURE_SUBSPACE(PSI,BITS) gives result of the measurement
%   in binary state description. The kth bit of BIN corresponds to the kth
%   measured qubit from the left. 
%
%   Note that measurement collapsed PSI to a random result PHI.
%   The probability distributions of other possible PHI's are not available
%   as outputs of this function.  This is intentional as this information
%   is lost in real life.  (However it is represented internally and the
%   user may peek at it in debug mode if interested.)

b = log2(size(psi,1));
assert(b==floor(b));

%convert to struct representation.
%each struct in the superposition array has a binary string state representation.
s=vec2struct(psi); %s is the struct representation

% change the order of the bits,such that
% bits to be measured are at the beginning
n = length(bits);          % number of qubits to measure
m = b - n;                 % number of qubits not measured
index = 1:b;
index(bits) = [];
index = [bits index];
index_inv(index) = 1:b; %this gives the inverse permutation to undo the sort
for iter=1:length(s)
    s(iter).bin = s(iter).bin(index);
end

%We can now follow the method in Gruska, p.71
psi = struct2vec(s);
l = 2^n;
k = 2^m;
pi = zeros(l,1);
for i=1:l
    pi(i) = sum((abs(psi((i-1)*k+(1:k)))).^2);
end

%simulate the collapse
[sp,ip]=sort(pi);
sp=cumsum(sp);     %the roulette wheel array
r=rand;            %pick a point on the roulette wheel
e=1;
while r>sp(e)      %find which part of wheel the point is in
    e=e+1;
end

%e is the index of the winning collapsed event
obs=ip(e);  
dec=obs-1;              % (adjust for matlab indices)
if nargout > 2
    bin=dec2bin(dec,n); % observed value of measured qubits
end
phi = zeros(k*l, 1);
phi((obs-1)*k+(1:k)) = (1/sqrt(pi(obs))).*psi((obs-1)*k+(1:k)); % collapsed state

%convert to struct representation.
s=vec2struct(phi); %s is the struct representation

%invert the sort, so the measured and unmeasured bits are back
%in their original places
for iter=1:length(s)
    s(iter).bin = s(iter).bin(index_inv);
end
phi = struct2vec(s);

