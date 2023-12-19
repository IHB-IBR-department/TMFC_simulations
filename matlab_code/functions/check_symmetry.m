function [r s] = check_symmetry(matrix)

% ========================================================================
% Ruslan Masharipov, October, 2023
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

matrix(1:1+size(matrix,1):end) = NaN;

up = matrix.';
up = up(tril(true(size(up)))).';  
up(isnan(up)) = [];
low = matrix;
low  = low(tril(true(size(low)))).';
low(isnan(low)) = [];

r = corr(low',up');
s = corr(low',up','Type','Spearman');

end