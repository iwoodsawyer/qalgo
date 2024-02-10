%% Solve system of linear equations
A = [1, -1/3; -1/3, 1]; 
b = [0; 1];
x = A\b;
x = x./norm(x);
prob = x.^2

% Get eigenvalues
n = 2;
t = 3*pi/(2^n);
[V,D] = eig(A);
d = diag(D);
U = V*diag(exp(d.*t.*1i))*V'; % = expm(A.*t.*1i);
[V,D] = eig(U);

v = d.*(2^n)*t/(2*pi); 
theta = 2*asin(1./v); 

% Find the expected ratio of measurement
probx = hll(U,b,n,theta)
