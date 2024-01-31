function Z = rotz(theta,n)
%ROTZ Create Rotation Z matrix
%  Z=ROTZ(THETA,N) Create N-bit Rotation Z matrix.

if n==1
    Z = [ 1 0; 0 exp(2*pi*1i*theta)];
else
    Z1 = rotz(theta,1);
    Z=1;
    for i=1:n
        Z=kron(Z,Z1);
    end
end

