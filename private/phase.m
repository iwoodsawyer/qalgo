function S = phase(n)
%Phase Create Phase S matrix
%   S=PHASE(N) Create N-bit PHASE matrix.

if n==1
    S = [ 1 0; 0 1i ];
else
    S1 = phase(1);
    S=1;
    for i=1:n
        S=kron(S,S1);    
    end
end

