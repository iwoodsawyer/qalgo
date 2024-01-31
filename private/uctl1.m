function [Uctl,psi] = uctl1(U, psi, k)
%UCTL Return the matrix representation of the Clt-U^(2^(k-1)) gate.

psi_u = U^(2^(k-1)) * psi;

for j = 1:size(psi, 1)
    if psi(j) ~= 0
        phase = psi_u(j) / psi(j);
    end
end

Uctl = [1 0; 0 phase];
