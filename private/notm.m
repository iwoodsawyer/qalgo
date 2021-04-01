function H = notm(n)
%NOTM Create NOT  matrix
%   H=NOTM(N) Create N-bit NOT matrix.

if n==1
    H = [ 0 1; 1 0 ];
else
    H1 = notm(1);
    H=1;
    for i=1:n
        H=kron(H,H1);    
    end
end

