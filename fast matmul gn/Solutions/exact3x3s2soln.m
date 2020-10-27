function [A,B,C,D] = exact3x3s2soln()

    % solution for 3x3 matmul with 2 symmetric components and 7 triplets
    U =[
         0  ,   1 ,    0 ,    0  ,   0  ,   0   ,  0   ,  0    , 0    , 0 ,   -1,     0   ,  0    , 0   ,  0  ,   1  ,   0   ,  0   ,  0   ,  0    , 1   , -1   ,  0;
         0  ,   0  ,   0  ,   1  ,   0   ,  0  ,   0  ,   0  ,   0     0  ,   0  ,   0  ,   0  ,   0  ,   0  ,  -1  ,   0  ,   0  ,   0  ,   1  ,  -1  ,   1   ,  0;
         0   ,  0  ,   0  ,   0  ,  -1   ,  0  ,   0  ,   0  ,  -1     0  ,   0  ,   0  ,  -1  ,   1  ,   0  ,   1  ,   0  ,   0  ,   0  ,  -1  ,   1  ,  -1    , 0;
         0   ,  0  ,   0  ,   0  ,   0  ,   0  ,   1  ,   0  ,   0     0  ,   0  ,   1  ,   1  ,   0  ,   0  ,   0  ,   0  ,  -1  ,   0  ,   0  ,   0  ,   0   ,  1;
         1  ,   0  ,   0  ,   1  ,   0  ,   0  ,   0  ,   0  ,  -1    -1  ,   0  ,   0  ,   0  ,   0  ,   1  ,   0  ,   1  ,   0  ,  -1  ,   0  ,   0  ,   0   ,  0;
         0  ,   0  ,   0  ,   0  ,   0  ,   0  ,  -1  ,   1  ,   0     1  ,   0  ,   0  ,  -1  ,   0  ,   0  ,   0  ,   0  ,   0  ,   0  ,   0  ,   0  ,   0   , -1;
         0  ,   0  ,   0  ,   0  ,   0  ,  -1  ,   1  ,   0  ,   0     0  ,  -1  ,   0  ,   0  ,   0  ,   0  ,   1  ,   0  ,   0  ,   0  ,   0  ,   0  ,  -1   ,  0;
         0  ,   0  ,   0  ,   0  ,   0  ,   0  ,   0  ,   0  ,   0     0  ,   0  ,   0  ,   0  ,   0  ,   0  ,   0  ,   1  ,   0  ,  -1  ,   0  ,   0  ,   1   ,  0;
         1  ,   0  ,   1  ,   0  ,   0  ,   0  ,  -1  ,   1  ,   0     0  ,   0  ,   0  ,   0  ,   0  ,   0  ,   0  ,   0  ,   0  ,   0  ,   0  ,   0  ,  -1   , -1;
    ];

    % separating U into ABCD matrices
    % note that C and D are swapped....
    A = U(:,1:2);
    B = U(:,3:9);
    D = U(:,10:16);
    C = U(:,17:23);

    % checking that solution is exact
    assert( norm(tensor_residual(A,B,C,D,matmul_tensor(3,3,3))) == 0 );

end