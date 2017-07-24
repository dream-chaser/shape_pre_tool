function [evecs,evals,area] = calc_lbo(shape, n)

[W, A] = calcLB(shape);
area = diag(A);

try
    [evecs,evals] = eigs(W, A, n, 'SM');
catch
    [evecs,evals] = eigs(W, A, n, -1e-5);
end
evals = abs(diag(evals));

% [evals, perm] = sort(evals);
% evecs = evecs(:, perm);
