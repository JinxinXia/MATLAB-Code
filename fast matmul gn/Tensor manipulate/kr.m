function AB = kr(A,B)
% conduct Khatri-Rao product
[m1,n1] = size(A);
[m2,n2] = size(B);

if n1 ~= n2
    error(' Error in kr - The matrices must have the same number of columns');
end

AB = zeros(m1*m2,n1);
   
for i = 1:n1
    AB(:,i) = kron(A(:,i),B(:,i));
end

end

