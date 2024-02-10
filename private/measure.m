function [phi,dec,bin] = measure(psi)
%MEASURE Measure phi wrt standard basis.
%   PHI=MEASURE(PSI) returns phi which is the new state of psi after
%   measurement collapse.
%
%   [PHI,DEC]=MEASURE(PSI) gives result in decimal state description.
%
%   [PHI,DEC,BIN]=MEASURE(PSI) gives result of the measurement in binary
%   state description.

b = log2(size(psi,1));
assert(b==floor(b));

p = norm(psi); 
assert(abs(sum(p)-1)<sqrt(eps), 'Psi not normalised.');

%simulate the collapse
[sp,ip]=sort(p);
sp=cumsum(sp);
r=rand;
i=1;
while r>sp(i)
    i=i+1;
end

%i is the index of the winning collapsed event
obs=ip(i);  
phi=zeros(size(psi));
phi(obs)=1;             % collapsed state
dec = obs-1;            % (adjust for matlab indices)
if nargout > 2
    b = log2(size(psi,1));
    bin = dec2bin(dec, b);  % observed value of measured qubits
end
