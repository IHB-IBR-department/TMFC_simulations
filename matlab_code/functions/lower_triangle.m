function low = lower_triangle(matrix)

% ========================================================================
% Ruslan Masharipov, October, 2023
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

matrix(1:1+size(matrix,1):end) = NaN;
low = matrix(tril(true(size(matrix)))).';
low(isnan(low)) = [];

end