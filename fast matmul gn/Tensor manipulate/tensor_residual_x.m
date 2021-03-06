function [r] = tensor_residual_x(x)
% x is the vectorization of matrices A, B, C, D

Tau = matmul_tensor(2,2,2);

if size(size(Tau),2) ~= 3
    error('The input must be a three way tensor');
end

% get the size for A and B, where B, C, D share the same column size
% in this function we set the S and T to get A,B,C,D 
I = size(Tau,1);
S = 1;
T = 2;


A = reshape(x(1:I*S),[S,I])';
B = reshape(x(I*S+1:I*S+I*T),[T,I])';
C = reshape(x(I*S+I*T+1:I*S+2*I*T),[T,I])';
D = reshape(x(I*S+2*I*T+1:I*S+3*I*T),[T,I])';

R = Tau;

% main for loop to create r 
for i = 1:I
    for j = 1:I
        for k = 1:I
            for s = 1:S
                R(i,j,k) = R(i,j,k) - A(i,s)*A(j,s)*A(k,s);
            end
            for t = 1:T
                R(i,j,k) = R(i,j,k) - B(i,t)*D(j,t)*C(k,t) - C(i,t)*B(j,t)*D(k,t) - D(i,t)*C(j,t)*B(k,t);
            end
        end
    end
end
                
r = R(:);             
               

end