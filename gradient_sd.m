% Gradient Steepest Descent Function
% --- Parameters ---
% A = N x N non-singular matrix
% b = N x 1 vector
% x0 = N x 1 vector
% --- Return Value ---
% xks = every xk value found along the way is stored in a list
function xks = gradient_sd(A, b, x0)

    % Initializing variables
    xks = [];
    xks = [xks; x0.'];
    rk = b - A*x0;
    xk = x0;
    k = 0;
    
    % Stepping through the algorithm
    while true
       
        % Calculating ak based on rk and A
        ak = (rk.'*rk) / (rk.'*A*rk);
        x = xk + ak*rk;
        k = k + 1;
        
        % Appending new xk to the list of xk values
        xks = [xks; x.'];
        
        % Calculating the error
        error = norm(x - xk);
        
        % If the algorithm has converged, exit
        if (error < 1.0e-7)
            return;
        end
        
        % Update variables
        rk = rk - ak*A*rk;
        xk = x;       
        
    end
end

