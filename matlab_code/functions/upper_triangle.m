function up = upper_triangle(matrix)

% ========================================================================
% Ruslan Masharipov, October, 2023
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

matrix(1:1+size(matrix,1):end) = NaN;
matrix = matrix.';
up = matrix(tril(true(size(matrix)))).';
up(isnan(up)) = [];

end