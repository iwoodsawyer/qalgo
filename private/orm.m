function H = orm(n)
%ORM Create OR matrix
%   H=ORM(N) Create N-by-2N-bit OR matrix.

if n==1
    H = [ 1 0 0 0; 0 1 1 1];
else
    H1 = orm(1);
    H=1;
    for i=1:n
        H=kron(H,H1);    
    end
end

