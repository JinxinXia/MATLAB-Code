% get exact solution for 3x3 case
[A,B,C,D] = exact3x3s2soln();

% perturb exact solution by relative eta
eta = 1e-4;
A = A + eta*randn(size(A)).*A;
B = B + eta*randn(size(B)).*B;
C = C + eta*randn(size(C)).*C;
D = D + eta*randn(size(D)).*D;
norm_before = norm(tensor_residual(A,B,C,D,tau))

% call GaussNewton to confirm convergence back to exact solution
tau = matmul_tensor(3,3,3);
tol = 1e-6;
max_step = 20;
lambda = 1e-1;
[A,B,C,D] = GaussNewton(tau,tol,max_step,lambda,A,B,C,D);

norm_r = norm(tensor_residual(A,B,C,D,tau))