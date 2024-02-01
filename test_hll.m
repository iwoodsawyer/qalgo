% Solve system of linear equations
A = [1, -1/3; -1/3, 1]; 
b = [0; 1];
x = A\b;
x = x./norm(x);
prob = x.^2

% Get eigenvalues
t = 3*pi/4;
[V,D] = eig(A);
d = diag(D);
U = V*diag(exp(d.*t.*1i))/V; % = expm(A.*t.*1i);
v = d/min(d); % scale such that 1/v <= 1

% Find the expected ratio of measurement
x = hll3(U,b,v)

