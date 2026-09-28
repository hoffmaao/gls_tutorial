function E = window_errors(p, W)
%WINDOW_ERRORS Reference solution, step 4: error bars and the chi-square check.
%
% E = ref.window_errors(p, W)
%
% The residuals are whitened with each row's full covariance, so the
% Jacobian J of r gives the posterior covariance directly (lessons 1, 3, 6):
%   C_M      = inv(J'J)
%   dof      = sum(W.rank) - 4        (independent numbers minus 3 + g)
%   chi2_dof = cost / dof
%   sigma    = sqrt(diag(C_M)) * sqrt(max(chi2_dof, 1))   widen, never shrink
%
% E fields: sigma [3 x 1], C_M [3 x 3], chi2_dof, dof, g

[cost, g, r] = ref.window_cost(p, W);
h = [1e-4; 1e-4; 1e-6];
J = zeros(numel(r), 3);
for j = 1:3
  pj = p; pj(j) = pj(j) + h(j);
  [~, ~, rj] = ref.window_cost(pj, W);
  J(:, j) = (rj - r) / h(j);
end
C_M = inv(J' * J);
dof = sum(W.rank) - 4;
chi2_dof = cost / dof;
E = struct('sigma', sqrt(max(diag(C_M), 0)) * sqrt(max(chi2_dof, 1)), ...
  'C_M', C_M, 'chi2_dof', chi2_dof, 'dof', dof, 'g', g);
end
