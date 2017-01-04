% Function to decompose rectangular matrix H into corresponding parts
% Q and R where Q is an orthoganal matrix and R is an upper triangular
% matrix
function [Q, R] = qr_decomp(H)

    % Initialization
    [m, n] = size(H);

    Q = zeros(m, n);
    R = zeros(n, n);
    V = H;

    % QR Decomposition Algorithm
    for j = 1:n
        
        R(j, j) = norm(V(:, j));
        Q(:, j) = (V(:, j) / R(j, j));
        
        for i = (j+1):n
           
            R(j, i) = dot(Q(:, j).', V(:, i));
            V(:, i) = V(:, i) - (R(j, i) * Q(:, j));
            
        end 
    end
end

