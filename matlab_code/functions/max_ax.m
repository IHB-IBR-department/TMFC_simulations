function min_max_ax = max_ax(matrix,type)

% ========================================================================
% Ruslan Masharipov, October, 2023
% email: ruslan.s.masharipov@gmail.com
% ========================================================================

matrix(1:1+size(matrix,1):end) = NaN;
max_ax = roundf(max(max(abs(matrix))));

if type == 1
    min_max_ax = [-max_ax max_ax];
elseif type == 0   
    min_max_ax = [-max_ax 0 max_ax];
end
    
end

function y = roundf(x)
y = round(x*10000)/10000;
end