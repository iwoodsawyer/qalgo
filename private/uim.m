function Ui = uim(f,n,y)
%UIM Make unitary f-conditional inverter UI matrix from function F
%   UI = UIM(F,N) makes unitary f-conditional inverter UF from function F,
%   where F must be of the form F(X,N) where X is a M-bits integer, and N
%   is the numbers of input and output bits.

l = 2^n;
Ui = speye(l);

if nargin < 3
    for i=1:l
        Ui(i,i) = (-1)^feval(f,i-1);
    end
else
    for i=1:l
        Ui(i,i) = (-1)^feval(f,i-1,y);
    end
end



