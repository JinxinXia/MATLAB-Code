function [A,B,C,D] = TLab_GaussNewton(Tau,tol,max_step,A,B,C,D)

X = [A(:); B(:); C(:); D(:)];

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


iter = 1;
while r_norm > tol && iter < max_step
    % compute residual
    r  = tensorResX(X);
    r_norm = norm(r);
    
    % set variable matrix
    X = [A(:);B(:);C(:);D(:)];
    
    % create target by imposing [-1,1] range on all variables
    target = X;
    bigpos = (X >  1);
    bigneg = (X < -1);
    target(bigpos) = 1;
    target(bigneg) = -1;
    
    % compute p using tensorlab Gauss Newton function
    [p,~] = nls_gndl(@tensorResX,@jacobian,X);
  
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





