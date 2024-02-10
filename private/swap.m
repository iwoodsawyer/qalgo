function S = swap(i,j,n)
%SWAP Swap qubits I and J in N qubits

I = identity(1);
SI = 1;
for k = n:-1:1
    if j == k
        SI = kron(SI,I);
    elseif i == k
        SI = kron(SI,I);
    else
        SI = kron(SI,I);
    end
end
SX = 1;
for k = n:-1:1
    if j == k
        SX = kron(SX,paulix(1));
    elseif i == k
        SX = kron(SX,paulix(1));
    else
        SX = kron(SX,I);
    end
end
SY = 1;
for k = n:-1:1
    if j == k
        SY = kron(SY,pauliy(1));
    elseif i == k
        SY = kron(SY,pauliy(1));
    else
        SY = kron(SY,I);
    end
end
SZ = 1;
for k = n:-1:1
    if j == k
        SZ = kron(SZ,pauliz(1));
    elseif i == k
        SZ = kron(SZ,pauliz(1));
    else
        SZ = kron(SZ,I);
    end
end

S = (SI + SX + SY + SZ)./2;