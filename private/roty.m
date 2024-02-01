function Y = roty(theta,n)
%ROTY Create Rotation Y matrix
%  Y=ROTY(THETA,N) Create N-bit Rotation Y matrix.

if n==1
    Y = [cos(theta/2),-sin(theta/2);...
         sin(theta/2),cos(theta/2)];
else
    Y1 = roty(theta,1);
    Y=1;
    for i=1:n
        Y=kron(Y,Y1);
    end
end

