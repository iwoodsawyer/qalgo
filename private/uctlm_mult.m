function [Im,Upow] = uctlm_mult(U, k, n, nb)
%UCTLM Return the matrix representation of the Clt-U^(2^(k-1)) gate.

I = identity(1);
I00 = I;
I00(2,2) = 0;
I11 = I;
I11(1,1) = 0;

Im = cell(n+1,1);
Im{1} = identity(nb);
for j = n:-1:1
    if j == k
        Im{n-j+2} = I00;
    else
        Im{n-j+2} = I;
    end
end

Upow = cell(n+1,1);
Upow{1} = U^(2^(k-1));
for j = n:-1:1
    if j == k
        Upow{n-j+2} = I11;
    else
        Upow{n-j+2} = I;
    end
end





