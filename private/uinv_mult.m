function [Im,Uinv] = uinv_mult(d, k, n)
%UINV Return the matrix representation of the Clt-Ry gate.

theta = 2*asin(1/d(k));
Ry = roty(theta,1);

I = identity(1);
I00 = I;
I00(2,2) = 0;
I11 = I;
I11(1,1) = 0;

Im = cell(n+1,1);
Im{n+1} = I;
for j = 1:n
    if j == k
        Im{n-j+1} = I00;
    else
        Im{n-j+1} = I;
    end
end

Uinv = cell(n+1,1);
Uinv{n+1} = Ry;
for j = 1:n
    if j == k
        Uinv{n-j+1} = I11;
    else
        Uinv{n-j+1} = I;
    end
end


