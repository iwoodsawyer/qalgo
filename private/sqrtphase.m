function T = sqrtphase(n)
%SQRTPHASE Create squareroot phase T matrix
%   T=SQRTPHASE(N) Create N-bit SQRTPHASE matrix.

if n==1
    T = [ 1 0; 0 exp((pi/4)*1i)];
else
    T1 = sqrtphase(1);
    T=1;
    for i=1:n
        T=kron(T,T1);    
    end
end

