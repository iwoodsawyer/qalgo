function [Uctl,psi] = uctl(U, psi, k, n)
%UCTL Return the matrix representation of the Clt-U^(2^(k-1)) gate.

psi_u = U^(2^(k-1)) * psi;

for j = 1:size(psi, 1)
    if psi(j) ~= 0
        phase = psi_u(j) / psi(j);
    end
end

% Uctl = 1;
% for j = n:-1:1
%     if j == k
%         Uctl = kron(Uctl, [1 0; 0 phase]);
%     else
%         Uctl = kron(Uctl, eye(2));
%     end
% end

Uctl = cell(n,1);
for j = 1:n
    if j == k
        Uctl{n-j+1} = [1 0; 0 phase];
    else
        Uctl{n-j+1} = eye(2);
    end
end

end