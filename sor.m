% SOR method for solving linear equations
% --- Parameters ---
% A = N x N non-singular matrix
% b = N x 1 vector
% x0 = N x 1 vector
% w = convergence factor
% --- Return Value ---
% xks = every xk value found along the way is stored in a list
function xks = sor(A, b, x0, w)

    % Initializing variables
    xks = [];
    xks = [xks; x0.'];
    xk = x0;
    xk_recent = x0;
    n = size(x0, 1);
    xk_new = zeros(n, 1);

    % Iteratively calculating a new x value as a solution to the
    % equation Ax = b
    while true
        
        % Calculating next x value
        for i = 1:n           
            temp = b(i);
            for j = 1:n            
                if (j ~= i)
                    temp = temp - A(i,j)*xk_recent(j);
                end                
            end

            xk_new(i) = (w/A(i,i))*temp + (1-w)*xk(i);
            
            % Updating most recent x value found since we just
            % calculated a new one
            xk_recent(i) = xk_new(i);
        
        end
        
        % Appendng current xk value to list
        xks = [xks; xk_new.'];
        
        % Calculating error
        error = norm(xk_new - xk);
        
        % If the algorithm has converged, exit
        if (error < 1.0e-7)
            return;
        end
        
        % Updating xk for the next step
        xk = xk_new;
        
    end
end

