
function J = jacobian(C, D)

     m = size(C,1);
     I = eye(m);

    
        % P2 and P3 are permutation matrices, P2 : ijk to jik ;   P3 : ijk to kij
        [r2,r3] = build_perm_matrices(m);
        
        K2 = kron(I,kr(C,D));
        K3 = kron(I,kr(D,C));
      
        
        J = kron(I,kr(D,C)) + K2(r2,:) + K3(r3,:);
   

end



