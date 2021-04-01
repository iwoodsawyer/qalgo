function D = iam(n)
%IAM Inversion about the average for N bits
%   D=IAM(N)

l = 2^n;
D = (2/l).*ones(l) - eye(l);