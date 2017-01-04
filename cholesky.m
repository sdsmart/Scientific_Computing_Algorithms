% This function performs the Cholesky algorithm for decomposing
% a given matrix A. A is assumed to be S.P.D.
%
% Returns a unit-lower triangular matrix G
% when multiplied by its transpose should compose
% back into A (with some small error).
function G = cholesky(A)

    % Calculating the dimensions of A
    % and A is a square matrix.
    n = size(A, 1);
    
    % Initializing G.
    G = zeros(n);
    
    % Beginning the Cholesky algorithm
    for j = 1:n
        
        % Computing the diagonal elements of G
        temp = 0;
        for k = 1:(j-1)
            temp = temp + G(j, k)^2;
        end
        G(j, j) = sqrt(A(j, j) - temp);
        
        % Recovering the columns of G
        for i = (j+1):n
            temp = 0;
            for k = 1:(j-1)
                temp = temp + G(i, k) * G(j, k);
            end
            G(i, j) = (A(i, j) - temp) / G(j, j);
        end 
    end
end

