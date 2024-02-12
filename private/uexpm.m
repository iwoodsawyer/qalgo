function [U,d,n,t] = uexpm(A)
%UEXPM Unitary Matrix Exponent
%[U,D,N,T] = UEXP(A) determines U = EXPM(A.*t.*1i) where T is calculated 
%  such that T<=(2*PI)*MAX(ABS(EIG(A))) and MIN(EIG(U)) is an integer. 

[V,D] = eig(A);
d = diag(D);
n = max(size(A,1),floor(log2(max(abs(d)))+1+sqrt(eps))); % check number of qubits to represent largest eigenvalue
d_min_approx = (0:(2^n - 1))'; 
t = (2*pi)/(2.^n).*d_min_approx./min(abs(d)); % calculate t that have interger
t = max(t(t < ((2*pi)/max(abs(d)) - sqrt(eps))));
U = V*diag(exp(d.*t.*1i))*V'; % = expm(A.*t.*1i);

