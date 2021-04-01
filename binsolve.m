function x = binsolve(A,b)
%BINSOLVE Solve system of binary equations
%  X = BINSOLVE(A,B) solves the system of binary equations AX=B.

mat = [A b];
n = size(A,1);
ma = size(A,2);
assert(n==ma);
mb = size(b,2);
m = ma+mb;

% Gauss Elimination with partial pivoting on the matrix
for i = 1:n
    if ~bin2dec(mat(i,i))
        for j = i+1:1:n;
            if bin2dec(mat(j,i))
                for k = i:m
                    t = mat(i,k);
                    mat(i,k) = mat(j,k);
                    mat(j,k) = t;
                end
            end
        end
    end
    for j = i+1:n
        if bin2dec(mat(j,i))
            for k = i:m
                mat(j,k) = dec2bin(bitxor(bin2dec(mat(j,k)),bin2dec(mat(i,k))));
            end
        end
    end
end

% Guass Jordan Elimination
for i = n:-1:1
    if bin2dec(mat(i,i))
        for j = i-1:-1:1
            if bin2dec(mat(j,i))
                for k = i:m
                    mat(j,k) = dec2bin(bitxor(bin2dec(mat(j,k)),bin2dec(mat(i,k))));
                end
            end
        end
    end
end

% Determine rank
nzi = find(bin2dec(mat(1:n,1:n)));
rank = length(nzi);

if rank == n
    x = mat(:,ma+1:m);
else
    error('Matrix A has low rank!');
end
