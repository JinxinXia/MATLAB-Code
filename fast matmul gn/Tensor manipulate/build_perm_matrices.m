function [r2,r3] = build_perm_matrices(n)

    % build vec permutation (swaps row for col major of matrix)
    VP = zeros(n^2);
    for i=1:n 
        for j=1:n
            H = zeros(n); 
            H(i,j)=1; 
            VP = VP + kron(H,H'); 
        end
    end
    
    % build P2 and P3
    P2 = kron(VP,eye(n));
    P3 = kron(eye(n),VP)*P2; 
   
    % using the row vector of transponsed permutation matrix is cheaper
    [r2,~,~] = find(P2');
    [r3,~,~] = find(P3');
    
   
    
end
