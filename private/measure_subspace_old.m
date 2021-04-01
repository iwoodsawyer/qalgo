function [phi,dec,bin] = measure_subspace(psi, bits)
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
if ~isstruct(psi) %if vec
    s=vec2struct(psi); %s is the struct representation
end

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
l = 2^n;
k = 2^m;
ic = zeros(length(s),1);
jc = zeros(length(s),1);
vc = zeros(length(s),1);
for iter=1:length(s)
    %adjust for 1-based matlab indices
    ic(iter) = bin2dec(s(iter).bin(1:n))+1;   %decimal rep of measured part
    jc(iter) = bin2dec(s(iter).bin(n+1:b))+1; %decimal rep of unmeasured part
    vc(iter) = s(iter).alpha; 
end
c = sparse(ic,jc,vc,l,k); %coefficients matrix
p_i = sum(abs(c.^2), 2)'; %probability distribution of observation (from 0 to n-1)

%We could pick a random outcome now, then just calculate its output state
%Instead we compute the output states for each outcome, as a matrix and
%choose the column corresponding to the randomly chosen state. This is more
%work, but it may be of interest to inspect the whole set of possible
%outcomes at this stage. (They are not returned as output variables, since
%in the real world his infomation is not available. The user can sneak a
%glipse at the mind of god by inspecting them in debug mode if required!)

phi_i = zeros(length(psi), l); %set of possible output states (concenated horizontally)
for i=1:l
    for j=1:k   %add up contributions from each j
        if (c(i,j)~=0)
            bin = [dec2bin(i-1,n) dec2bin(j-1,m)];
            
            %invert the sort, so the measured and unmeasured bits are back
            %in their original places
            bin_inv = bin(index_inv);
            
            phi_i(:,i) = phi_i(:,i) + c(i,j) * bin2vec(bin_inv);
        end
    end
    if (p_i(i)~=0)
        phi_i(:,i) = phi_i(:,i) / sqrt(p_i(i)); %normalise
    end
end

%simulate the collapse
[p_sort,p_index]=sort(p_i);
sp=cumsum(p_sort);          %the roulette wheel array
r=rand;                     %pick a point on the roulette wheel
e=1;
while r>sp(e)               %find which part of wheel the point is in
    e=e+1;
end

%e is the index of the winning collapsed event
obs = p_index(e);
phi = phi_i(:,obs);     % collapsed state
dec = obs-1;            % adjust for matlab indices
if nargout > 2
    bin = dec2bin(dec, n);  % observed value of measured qubits
end