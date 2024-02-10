function Uctl = uctl1(U, k, nb)
%UCTL1 Return the matrix representation of the Clt-U^(2^(k-1)) gate.

I = identity(1);
I11 = I;
I11(1,1) = 0;
Im = identity(nb);
Uctl = U^(2^(k-1)) - Im;
Uctl = kron(Uctl,I11);


