% here I should put 2 2 2 for N K M 
function T = matmul_tensor(N,K,M) 
% build N x K x M matmul tensor
% assumes col-major order of inputs and row-major order of output
 

    % init tensor
    %T = tenzeros(N*K,K*M,N*M);
    

    % iterate over nonzero entries in tensor
    for k = 1:K,
        for m = 1:M,
            for n = 1:N,
                T(N*(k-1)+n,K*(m-1)+k,M*(n-1)+m) = 1;
            end
        end
    end
    

end
