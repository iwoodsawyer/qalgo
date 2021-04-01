function [rank,mat] = binrank(mat)
%BINRANK Determine rank of binary matrix
%  RANK = BINRANK(MAT) determines rank of binary matrix MAT.

n = size(mat,1);
m = size(mat,2);

% Gauss Elimination with partial pivoting on the matrix
for i = 1:min(m,n)
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
for i = min(m,n):-1:1
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
nzi = find(bin2dec(mat));
rank = length(nzi);

% return only independent rows
mat = mat(nzi,:);

