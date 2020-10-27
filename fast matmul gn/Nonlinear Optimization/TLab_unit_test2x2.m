
% get exact solution for 2x2 case
[A,B,C,D] = exact2x2s1soln();

% perturb exact solution by relative eta
eta = 1e-1;
A = A + eta*randn(size(A)).*A;
B = B + eta*randn(size(B)).*B;
C = C + eta*randn(size(C)).*C;
D = D + eta*randn(size(D)).*D;

% call GaussNewton to confirm convergence back to exact solution
tau = matmul_tensor(2,2,2);
tol = 1e-6;
max_step = 200;
lambda = 1e-2;
[A,B,C,D] = TLab_GaussNewton(tau,tol,max_step,A,B,C,D);

norm_r = norm(tensor_residual(A,B,C,D,tau))




