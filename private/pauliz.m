function Z = pauliz(n)
%PAULIZ Create Pauli Z matrix
%   Z=PAULIZ(N) Create N-bit PAULIZ matrix.

if n==1
    Z = [ 1 0; 0 -1 ];
else
    Z1 = pauliz(1);
    Z=1;
    for i=1:n
        Z=kron(Z,Z1);    
    end
end

