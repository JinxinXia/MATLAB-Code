function [A,B,C,D] = LevenbergMarquardt(Tau,tol,max_step,A,B,C,D)

if size(size(Tau),2) ~= 3
    error('First input must be a three way tensor');
end


% determine dimensions from inputs
I = size(Tau,1);
S = 1;
T = 2;

% determine initial residual
r = tensor_residual(A,B,C,D,Tau);
r_norm = norm(r);




% lsqcurvefit is the built in function for Levenberg Marquardt algorithm

iter = 1;

while r_norm > tol && iter < max_step
    
    % set variable matrix
    X = [A(:);B(:);C(:);D(:)];
    
    X_size = size(X,1);
    
    Xdata = zeros(X_size,10);
    Ydata = zeros(X_size,1);
    for i = 1:10
        Xdata(:,i) = X + 0.1*i;
        Ydata(:,i) = norm(tensor_residual_x(Xdata(:,i)));
    end
    
    % using Levenberg Marquardt
    options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
    lb = [];
    ub = [];
    p = lsqcurvefit('tensor_residual_x',X, Xdata, Ydata,lb,ub,options) - X;
    
  
    % extract updates to ABCD
    A = A + reshape(p(1:I*S),[S,I])';
    B = B + reshape(p(I*S+1:I*S+I*T),[T,I])';
    C = C + reshape(p(I*S+I*T+1:I*S+2*I*T),[T,I])';
    D = D + reshape(p(I*S+2*I*T+1:I*S+3*I*T),[T,I])';
        
    X = p + X;
    
    iter = iter + 1;
    r_norm = tensor_residual_x(X);
    
    if mod(iter,1) == 0
        fprintf("%1.3g, %1.3g\n",r_norm, norm(p));
    end
    
end % end main while loop



end