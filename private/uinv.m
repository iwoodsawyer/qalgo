function Uinv = uinv(d, k, n)
%UINV Return the matrix representation of the Clt-Ry gate.

theta = 2*asin(1/d(k));
Ry = roty(theta,1);

I00 = diag([1 0]);
I11 = diag([0 1]);
I = eye(2);

Im = eye(2);
for j = 1:n
    if j == k
        Im = kron(I00,Im);
    else
        Im = kron(I,Im);
    end
end

Uinv = Ry;
for j = 1:n
    if j == k
        Uinv = kron(I11,Uinv);
    else
        Uinv = kron(I,Uinv);
    end
end

Uinv = Im + Uinv;
