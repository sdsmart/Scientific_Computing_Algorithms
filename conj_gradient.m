% Conjugate Gradient Algorithm
% --- Parameters ---
% A = N x N non-singular matrix
% b = N x 1 vector
% x0 = N x 1 vector
% --- Return Value ---
% xks = every xk value found along the way is stored in a list
function xks = conj_gradient(A, b, x0)

    % Initializing variables
    xks = [];
    xks = [xks; x0.'];
    rk = b - A*x0;
    pk = rk;
    xk = x0;
    n = size(x0, 1);
    
    % Stepping through the algorithm
    for k = 0:(n-1)
       
        % Calculating ak based on rk, A, and pk
        ak = (rk.'*rk) / (pk.'*A*pk);
        x = xk + ak*pk;
        
        % Appending xk value to the list of xks
        xks = [xks; x.'];
        
        % Calculating rk and maintaining previous rk
        rk_old = rk;
        rk = rk - ak*A*pk;
        
        % Calculating the error
        error = rk.'*rk;
        
        % If the algorithm has converged, exit
        if (error < 1.0e-7)
            return;
        end
        
        % Calculating the new conjugate direction
        bk = (rk.'*rk) / (rk_old.'*rk_old);
        pk = rk + bk*pk;
        
        % Updating xk
        xk = x;
        
    end
end
