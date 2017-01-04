% This function performs the LU-Decomposition algorithm
% of a given matrix A. A is assumed to be non-singular.
%
% Returns a unit-lower triangular matrix L
% and an upper triangular matrix U.
function [L, U] = lu_decomp(A)

    % Calculating the dimensions of A
    % and A is a square matrix.
    n = size(A, 1);
    
    % Initializing L and U.
    L = zeros(n) + eye(n);
    U = zeros(n);
    
    % Beginning the LU-Decomposition Algorithm.
    for r = 1:n
        
        % Recovering the rows of U.
        for i = r:n
            temp = 0;
            for j = 1:(r-1)
                temp = temp + L(r, j) * U(j, i);
            end
            U(r, i) = A(r, i) - temp;
        end
        
        % Recovering the columns of L.
        for i = (r+1):n
            temp = 0;
            for j = 1:(r-1)
                temp = temp + L(i, j) * U(j, r);
            end
            L(i, r) = (A(i, r) - temp) / U(r, r);
        end        
    end    
end

