function Uinv = uinv(theta, k, n)
%UINV Return the matrix representation of the Clt-Ry gate.

Ry = roty(theta,1);
I = identity(1);
I11 = I;
I11(1,1) = 0;
Uinv = Ry - I;
for j = 1:n
    if j == k
        Uinv = kron(I11,Uinv);
    else
        Uinv = kron(I,Uinv);
    end
end
