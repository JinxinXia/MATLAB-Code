function [num_converg,norm_r,A,B,C,D]=fastmtm(N,randomseed,eta,max_step,tol,lambda)
%this function aims to find the fast matrix multiplication of N by N
%matrices, 


rng(randomseed);


Tau = matmul_tensor(N,N,N) ;
%[A,B,C,D,U,V,W] = strassen();


%I = size(Tau,1);

%choice of S = 1, T = 2 or S = 2, T = 7


%A = randn(I,S);
%B = randn(I,T);
%C = randn(I,T);
%D = randn(I,T);
%eta = 1e-2;
%A =  (eta*randn(size(A)));
%B =  (eta*randn(size(B)));
%C =  (eta*randn(size(C)));
%D = (eta*randn(size(D)));

 
%also need to check the change of the residuals

%max_step = 200;
%tol = 1e-5;
%lambda = 0.8;
%[solu] = GaussNewton(Tau,tol,max_step,lambda,A,B,C,D);

%rng(s);







num_converg = 0;
for i = 1:1
    
    %fprintf('\n%3g ',i);
    


    
A =  A - (eta*A);
B = B - (eta*B);
C =  C -  (eta*C);
D = D - (eta*D);

    
A =  (eta*randn(size(A)));
B =  (eta*randn(size(B)));
C =  (eta*randn(size(C)));
D =  (eta*randn(size(D)));

[A,B,C,D,iExit] = GaussNewton(Tau,tol,max_step,lambda,A,B,C,D);


num_converg = num_converg + iExit;

norm_r = norm(tensor_residual(round(A),round(B),round(C),round(D),Tau));
end


end