% This function performs the upper triangular
% system computation algorithm for a given matrix
% U and vector g. U should be an upper triangular matrix.
%
% Returns a solution vector x solving the equation
% Ux = g where U is a upper triangular matrix.
function x = upper_tri(U, g)

    % Calculating the number of rows of g.
    n = size(g, 1);
    
    % Initializing x.
    x = zeros(n, 1);
    
    % Beginning the upper triangular system
    % computation algorithm.
    x(n) = g(n) / U(n, n);
    for i = (n-1):-1:1
        temp = 0;
        for j = (i+1):n
            temp = temp + U(i, j) * x(j);
        end
        x(i) = (g(i) - temp) / U(i, i);
    end
end

