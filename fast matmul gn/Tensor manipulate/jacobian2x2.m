function J = jacobian2x2(X)
% return the jacobian with vectorized A, B, C, D as input the default I, S, 
% and T are the size of A, B, C, D. I is the row size and S is the column 
% size of A, and T is the column size of B, C, D
% poential size of I, S, T
% I = 4 , S = 1, T = 2
    I = 4;
    S = 1;
    T = 2;

    A = X(1:I*S);
    B = X(I*S+1:I*T+I*S);
    C = X(I*T+I*S+1:I*T+2*I*S);
    D = X(I*T+2*I*S+1:end);
    
    A = reshape(A,[I,S]);
    B = reshape(B,[I,T]);
    C = reshape(C,[I,T]);
    D = reshape(D,[I,T]);
    
    Ja = jacobian(A,A);
    Jb = jacobian(D,C);
    Jc = jacobian(B,D);
    Jd = jacobian(C,B);
    J = [Ja Jb Jc Jd];



end