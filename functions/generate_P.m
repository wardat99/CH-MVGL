function [P] = generate_P(n)
one_n = ones(n, 1); % Vector of ones

% Generate a random n x (n-1) matrix
A = rand(n, n-1);

% Apply Gram-Schmidt process to A to ensure orthogonality to the vector of ones
for i = 1:(n-1)
    % Subtract projection of A(:, i) onto one_n from A(:, i)
    A(:, i) = A(:, i) - (dot(A(:, i), one_n) / dot(one_n, one_n)) * one_n;
    
    % Orthogonalize against previous vectors in A
    for j = 1:(i-1)
        A(:, i) = A(:, i) - (dot(A(:, i), A(:, j)) / dot(A(:, j), A(:, j))) * A(:, j);
    end
    
    % Normalize the vectors to form an orthonormal basis
    A(:, i) = A(:, i) / norm(A(:, i));
end

% The resulting matrix A is your matrix P
P = A';
end