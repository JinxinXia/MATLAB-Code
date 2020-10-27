function [A,B,C,D] = GaussNewton(Tau,tol,max_step,lambda,A,B,C,D)


if size(size(Tau),2) ~= 3
    error('First input must be a three way tensor');
end

% determine dimensions from inputs
I = size(Tau,1);
S = size(A,2);
T = size(B,2);

% determine initial residual
r = tensor_residual(A,B,C,D,Tau);
r_norm = norm(r);

% initialize search direction
p = zeros(I*(3*T+S),1);

iter = 1;
while r_norm > tol && iter < max_step
    % build Jacobian 
    Ja = jacobian(A,A);
    Jb = jacobian(D,C);
    Jc = jacobian(B,D);
    Jd = jacobian(C,B);
    J = [Ja Jb Jc Jd];
    
    % compute residual
    r  = tensor_residual(A,B,C,D,Tau);
    r_norm = norm(r);
    
    % set variable matrix
    X = [A(:);B(:);C(:);D(:)];
    
    % create target by imposing [-1,1] range on all variables
    target = X;
    bigpos = (X >  1);
    bigneg = (X < -1);
    target(bigpos) = 1;
    target(bigneg) = -1;
    
    % solve the normal eqs for min_p ||Jp-r|| + lambda||p-(target-X)||
    % x_{k+1} = x_k + p, and we want 
    %                    x_{k+1}-target = p-(target-x_k) to be small
    coeff_matrix = J'*J + lambda*eye(length(p));
    rhs = J'*r + lambda*(target - X);
    p = coeff_matrix \ rhs;
    
    % extract updates to ABCD
    A = A + reshape(p(1:I*S),[S,I])';
    B = B + reshape(p(I*S+1:I*S+I*T),[T,I])';
    C = C + reshape(p(I*S+I*T+1:I*S+2*I*T),[T,I])';
    D = D + reshape(p(I*S+2*I*T+1:I*S+3*I*T),[T,I])';
        
    % balance largest entries of B, C, D
    mB = max(abs(B));
    mC = max(abs(C));
    mD = max(abs(D));
    
    B = B .* (ones(I,1)*(mC.*mD./(mB.^2)).^(1/3));
    C = C .* (ones(I,1)*(mB.*mD./(mC.^2)).^(1/3));
    D = D .* (ones(I,1)*(mB.*mC./(mD.^2)).^(1/3));
    
    iter = iter + 1;
    
    if mod(iter,1) == 0
        fprintf("%1.3g, %1.3g\n",r_norm, norm(p - (target-X)));
    end
    
end % end main while loop





end





