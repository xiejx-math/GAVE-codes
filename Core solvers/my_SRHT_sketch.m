function [SA, SB, Sb] = my_SRHT_sketch(A, B, b, ell)
    % Subsampled Randomized Hadamard Transform (SRHT) Sketch
    % Input:
    %   A: m x n matrix
    %   B: m x n matrix
    %   b: m x 1 vector
    %   ell: target dimension (ell < m)
    % Output:
    %   SA: ell x n compressed matrix
    %   SB: ell x n compressed matrix
    %   Sb: ell x 1 compressed vector
    % Coded by Jiaxin Xie, Beihang University, xiejx@buaa.edu.cn

    % Get the size of matrix A
    [m, n] = size(A);  

    % Check if ell is less than m
    if ell >= m
        error('Target dimension ell must be less than m.');
    end

    % Step 1: Randomly permute the rows of matrix A, B, and vector b
    % Generate a random permutation of row indices
    perm = randperm(m);  
    % Permute the rows of A, B, and b
    A_perm = A(perm, :);  
    B_perm = B(perm, :);  
    b_perm = b(perm);  

    % Step 2: Compute the Hadamard transform
    % Find the smallest power of 2 greater than or equal to m
    M = 2^nextpow2(m);  
    % Generate the Hadamard matrix of size M x M
    H = hadamard(M);  
    % Truncate the Hadamard matrix to the first m rows
    H = H(1:m, :);  

    % Step 3: Create a random sign matrix
    % Generate a random sign vector with entries +1 or -1
    D = randi([0, 1], m, 1) * 2 - 1;  
    % Convert the sign vector into a diagonal matrix
    D = diag(D);  

    % Step 4: Compute the SRHT for A, B, and b
    % Compute the Hadamard transform of the permuted and signed matrices
    SA = (1 / sqrt(M)) * H * D * A_perm;  
    SB = (1 / sqrt(M)) * H * D * B_perm;  
    Sb = (1 / sqrt(M)) * H * D * b_perm;  

    % Step 5: Randomly sample ell rows
    % Randomly select ell row indices
    sample_indices = randperm(m, ell);  
    % Sample the rows to get the compressed matrices and vector
    SA = SA(sample_indices, :);  
    SB = SB(sample_indices, :);  
    Sb = Sb(sample_indices);  

    % Step 6: Adjust the scaling factor
    % Scale the matrices and vector to preserve the energy
    SA = SA * sqrt(m / ell);  
    SB = SB * sqrt(m / ell);  
    Sb = Sb * sqrt(m / ell);  
end
