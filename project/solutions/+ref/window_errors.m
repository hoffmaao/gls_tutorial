function E = window_errors(p, W)
%WINDOW_ERRORS Reference solution, step 4: error bars and the chi-square check.
%
% E = ref.window_errors(p, W)
%
% Linearize at the answer (lesson 6, section 7): with J the Jacobian of the
% whitened residual, the posterior covariance is C_M = inv(J'J) (lesson 1,
% with the whitening of lesson 2 already inside J). Then (lesson 4):
%   dof      = number of real residuals - 3 nonlinear - 1 for g
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
dof = numel(r) - 4;
chi2_dof = cost / dof;
E = struct('sigma', sqrt(max(diag(C_M), 0)) * sqrt(max(chi2_dof, 1)), ...
  'C_M', C_M, 'chi2_dof', chi2_dof, 'dof', dof, 'g', g);
end
