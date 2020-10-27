function [r] = tensor_residual(A,B,C,D,Tau)


if size(size(Tau),2) ~= 3
    error('The input must be a three way tensor');
end


I = size(Tau,1);


% get the size for A and B, where B, C, D share the same column size
S = size(A,2);
T = size(B,2);

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