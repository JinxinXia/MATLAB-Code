function [A,B,C,D] = exact2x2s1soln() 
% return Strassen's ABCD representation

    % U
    U = zeros(4,7);
    U(1,1)=1; U(1,3)=1; U(1,5)=1; U(1,6)=-1;
    U(2,2)=1; U(2,6)=1;
    U(3,5)=1; U(3,7)=1;   
    U(4,1)=1; U(4,2)=1; U(4,4)=1; U(4,7)=-1;
    
    % ABCD
    A = U(:,1);
    B = [U(:,2) U(:,4)];
    C = [U(:,7) U(:,5)];
    D = [U(:,3) U(:,6)];
    
    % check solution
    assert( norm(tensor_residual(A,B,C,D,matmul_tensor(2,2,2))) == 0 );
    

end
