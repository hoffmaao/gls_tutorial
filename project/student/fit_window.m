function [p, cost, info] = fit_window(W, opts)
%FIT_WINDOW Step 3: the best [theta0; delta0; ddelta] for one window.
%
% [p, cost, info] = fit_window(W, opts)
%
% A. GRID SEARCH over theta0 (0:5:175 deg), ddelta (0 to 0.08 rad/m) and
%    delta0 (0:30:330 deg); keep the smallest window_cost. The misfit has
%    false minima (lesson 6), which a grid does not fall into.
%
% B. GAUSS-NEWTON from the best grid point (lesson 6):
%      r  = residuals at p (window_cost's third output)
%      J  = Jacobian: nudge each parameter by h(j), J(:,j) = (r_nudged - r)/h(j)
%      dp = -(J'*J) \ (J'*r)
%    If p + dp does not lower the cost, halve the step (up to 20 times).
%    Stop when no step helps or the cost changes by less than 1e-9 of itself.
%
% C. CONVENTION. (theta0 + 90 deg, -delta0, -ddelta) gives the same field.
%    Report the one with ddelta >= 0, then wrap theta0 into [0, pi) with
%    mod and delta0 into (-pi, pi] with angle(exp(1i*delta0)).
%
% opts (optional): .theta_step_deg (5), .ddelta_max (0.08), .n_ddelta (17),
%   .n_delta0 (12), .max_iter (30), .p0 (skip the grid, start here).
%
% Check:  check_step(3)

if nargin < 2, opts = struct(); end
th_g = deg2rad(0:opt(opts, 'theta_step_deg', 5):179);
dd_g = linspace(0, opt(opts, 'ddelta_max', 0.08), opt(opts, 'n_ddelta', 17));
nd0  = opt(opts, 'n_delta0', 12);
d0_g = (0:nd0-1) * 2*pi / nd0;

% --- A. grid search
if isfield(opts, 'p0') && ~isempty(opts.p0)
  p = opts.p0(:);
  cost = window_cost(p, W);
else
  cost = inf; p = [0; 0; 0];
  for th = th_g
    for dd = dd_g
      for d0 = d0_g
        % TODO 3a: evaluate this grid point; keep it if it beats the best so far.
        error('project:todo', 'fit_window grid search is not written yet - see TODO 3a');
      end
    end
  end
end
info.grid_cost = cost;

% --- B. Gauss-Newton with a halving line search
h = [1e-4; 1e-4; 1e-6];                     % nudges for theta0, delta0, ddelta
for it = 1:opt(opts, 'max_iter', 30)
  % TODO 3b: r, J, dp.
  error('project:todo', 'fit_window Gauss-Newton is not written yet - see TODO 3b');

  % TODO 3c: line search and stopping rules.
end
info.n_iter = it;

% --- C. convention
% TODO 3d: if p(3) < 0, switch to [p(1) + pi/2; -p(2); -p(3)]; then wrap.
end

function v = opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
