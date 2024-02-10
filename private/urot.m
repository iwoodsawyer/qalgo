function Urot = urot(theta, k, n)
%UROT Return the matrix representation of the Clt-Rz gate.

Rz = rotz(theta,1);
I = identity(1);
I11 = I;
I11(2,2) = 0;
Urot = Rz - I;
for j = 1:n
    if j == k
        Urot = kron(Urot,I11);
    else
        Urot = kron(Urot,I);
    end
end
