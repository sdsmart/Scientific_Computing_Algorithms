% This function performs the lower triangular
% system computation algorithm for a given matrix
% L and vector b. L should be a lower triangular matrix.
%
% Returns a solution vector g solving the equation
% Lg = b where L is a lower triangular matrix.
function g = lower_tri(L, b)

    % Calculating the number of rows of b.
    n = size(b, 1);
    
    % Initializing g.
    g = zeros(n, 1);
    
    % Beginning the lower triangular system
    % computation algorithm.
    g(1) = b(1) / L(1, 1);
    for i = 2:n
        temp = 0;
        for j = 1:(i-1)
            temp = temp + L(i, j) * g(j);
        end
        g(i) = (b(i) - temp) / L(i, i);
    end
end