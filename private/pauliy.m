function Y = pauliy(n)
%PAULIY Create Pauli Y matrix
%   Y=PAULIY(N) Create N-bit PAULIY matrix.

if n==1
    Y = [ 0 -1i; 1i 0 ];
else
    Y1 = pauliy(1);
    Y=1;
    for i=1:n
        Y=kron(Y,Y1);    
    end
end

