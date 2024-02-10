function Uctlm = uctlm(U, k, n, nb)
%UCTLM Return the matrix representation of the Clt-U^(2^(k-1)) gate.

I = identity(1);
I11 = I;
I11(1,1) = 0;
Im = identity(nb);
Uctlm = U^(2^(k-1)) - Im;
for j = n:-1:1
    if j == k
        Uctlm = kron(Uctlm,I11);
    else
        Uctlm = kron(Uctlm,I);
    end
end




