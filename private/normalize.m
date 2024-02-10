function psi = normalize(psi)
%NORMALIZE Normalize such that sum(phi.^2)=1
%   PSI=NORMALIZE(PSI)

psi = psi./norm(psi);