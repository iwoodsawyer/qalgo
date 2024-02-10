function Uctlm = uctlm_mult(U, k, n, nb)
%UCTLM Return the matrix representation of the Clt-U^(2^(k-1)) gate.

I = identity(1);
I11 = I;
I11(1,1) = 0;
Im = identity(nb);
Uctlm = cell(n+1,1);
Uctlm{1} = U^(2^(k-1)) - Im;
for j = n:-1:1
    if j == k
        Uctlm{n-j+2} = I11;
    else
        Uctlm{n-j+2} = I;
    end
end





