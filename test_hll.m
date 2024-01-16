A = [1, -1/3; -1/3, 1] 
b = [0; 1]

% Inverting the Matrix to find x %
inv_A = inv(A)
x = inv_A*b

% Find the expected ratio of measurement
norm_factor = norm(x)
x_norm = x/norm_factor
Prob_0 = x_norm(1,1)^2
Prob_1 = x_norm(2,1)^2

% Initiazing t and finding the eigenvalues and vector of the Matrix
t = 3*pi/4
[v,e] = eig(A)
V_tran = transpose(v)
Md = v*A*V_tran

% Diaginalizing the vectors needed for QPE rotation
y_diagonal = exp((i*Md*t));
y_diagonal(1,2) = 0;
y_diagonal(2,1) = 0;
y_non_diagonal = V_tran*y_diagonal*v

y_diagonal_2 = exp((i*Md*t)*2);
y_diagonal_2(1,2) = 0;
y_diagonal_2(2,1) = 0;
y_non_diagonal_2 = v*y_diagonal_2*V_tran

hll(A,b)
