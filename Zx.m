function y = Zx(W, x, y_ini)
%
% Estimates y = Z*x, with Z = inv(I - W)
% Use iterative method: the vector y must be such that Wy + x = y
%
% INPUT : 
% W is a n*n square matrix
% x is a n*1 column matrix
% y_ini is n*1 column matrix used as initial point for the iterative
% method.
%
% OUTPUT :
% y is a n*1 column matrix.
%
% if number of iterations exceeds MAX_ITER, an error is raised
%
% Robin Devooght : 2013, october 4th

PRECISION = realmin;
MAX_ITER = 1e3;

[n, m] = size(W);
[nx, mx] = size(x);
[ny, my] = size(y_ini);

if n ~= m
    error('Zx:W_square', 'W must be a square matrix');
end
if n~=nx
    error('Zx:Wx_dim', 'W and x dimensions must agree');
end
if n~=ny
    error('Zx:Wy_dim', 'W and y_ini dimensions must agree');
end
if mx~=1
    error('Zx:x_dim', 'x must be a column vector');
end
if my~=1
    error('Zx:y_dim', 'y_ini must be a column vector');
end

next_y = W*y_ini+x;
y = zeros(n,1);

for i = 1:MAX_ITER
    y = next_y;
    next_y = W*y+x;
    if max(abs(next_y - y)) < PRECISION
        % when precision is reached, the loop ends.
        break;
    end
end

if i == MAX_ITER 
    % if the maximum number of iterations has been reached, an error is raised.
	error('Zx:iteration_exceed', 'number of iterations reached MAX_ITER');
end

y = next_y;

if max(isinf(y)) == 1
    % raise an error if infinite value is found
    error('Zx:out_of_range', 'Iteration produced infinite values');
end
if max(isnan(y)) == 1
    % raise an error if NaN is found
    error('Zx:nan', 'Iteration produced NaN''s');
end

