function X = paulix(n)
%PAULIX Create Pauli X matrix
%   X=PAULIX(N) Create N-bit PAULIX matrix.

if n==1
    X = [ 0 1; 1 0 ];
else
    X1 = paulix(1);
    X=1;
    for i=1:n
        X=kron(X,X1);    
    end
end

