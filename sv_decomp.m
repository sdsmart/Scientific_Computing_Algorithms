% Function to decompose rectangular matrix H into corresponding parts
% U, E, and V where U is an mxn matrix, E is the square root of the eigen
% values of the grammian H.'*H and V is the transposed matrix of the eigen
% vectors of the grammian H.'*H
function [U, E, V] = sv_decomp(H)

    % Initialization
    [m, n] = size(H);
    
    U = zeros(m, n);

    % Utilizing MATLAB to find eigen values and vectors easily
    [V, E] = eig(H.'*H);
    E = sqrt(E);
    
    % Constructing the U matrix
    for i = 1:n
       
        U(:, i) = (1 / E(i, i)) * H * V(:, i);
        
    end
    
    % Transposing V for easy testing after values are returned
    V = V.';
end

