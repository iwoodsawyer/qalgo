function Uinv = uinv_mult(theta, k, n)
%UINV Return the matrix representation of the Clt-Ry gate.

Ry = roty(theta,1);
I = identity(1);
I11 = I;
I11(1,1) = 0;
Uinv = cell(n+1,1);
Uinv{n+1} = Ry - I;
for j = 1:n
    if j == k
        Uinv{n-j+1} = I11;
    else
        Uinv{n-j+1} = I;
    end
end


