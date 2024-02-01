function Uctlm = uctlm(U, k, n, nb)
%UCTLM Return the matrix representation of the Clt-U^(2^(k-1)) gate.

I00 = diag([1 0]);
I11 = diag([0 1]);
I = eye(2);

Im = eye(2^nb);
for j = n:-1:1
    if j == k
        Im = kron(Im,I00);
    else
        Im = kron(Im,I);
    end
end

Upow = U^(2^(k-1));
for j = n:-1:1
    if j == k
        Upow = kron(Upow,I11);
    else
        Upow = kron(Upow,I);
    end
end

Uctlm = Im + Upow;



