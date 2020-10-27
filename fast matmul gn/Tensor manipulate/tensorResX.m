function r = tensorResX(X)
% return tensor residuals that wiht X, vectorized A,B,C,D, as input

    % select size of A,B,C,D
    I = 4;
    S = 1;
    T = 2;

    % convert back to matrices
    A = X(1:I*S);
    B = X(I*S+1:I*T+I*S);
    C = X(I*T+I*S+1:I*T+2*I*S);
    D = X(I*T+2*I*S+1:end);
    
    % get the tensor for resiudal calculation
    Tau = matmul_tensor(2,2,2);
    
    % r is a residual vector
    r = tensor_residual(A,B,C,D,Tau);


end