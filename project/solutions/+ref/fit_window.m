function [p, cost, info] = fit_window(W, opts)
%FIT_WINDOW Reference solution, step 3: find the best [theta0; delta0; ddelta].
%
% [p, cost, info] = ref.fit_window(W, opts)
%
% 1. GRID SEARCH over all three (lesson 1, section 2 and lesson 6): cheap,
%    and it cannot be fooled by the false bottoms of a periodic model.
% 2. GAUSS-NEWTON from the best grid node (lesson 6): Jacobian of the
%    whitened residual by finite differences, step dp = -(J'J) \ (J'r),
%    halved until the cost drops.
% 3. CONVENTION. (theta0 + 90, -delta0, -ddelta) gives the same field, so
%    ddelta >= 0 is enforced: theta0 is the axis along which the phase
%    difference GROWS with depth (ptt.quadpolFabricLS does the same).
%
% opts: .theta_step_deg (5), .ddelta_max (0.08 rad/m), .n_ddelta (17),
%       .n_delta0 (12), .max_iter (30), .p0 (start here, skip the grid)

if nargin < 2, opts = struct(); end
th_g = deg2rad(0:opt(opts, 'theta_step_deg', 5):179);
dd_g = linspace(0, opt(opts, 'ddelta_max', 0.08), opt(opts, 'n_ddelta', 17));
nd0  = opt(opts, 'n_delta0', 12);
d0_g = (0:nd0-1) * 2*pi / nd0;

% --- 1. grid, unless a starting point was given (the jackknife, step 5,
% restarts from the full-data answer instead of searching again)
if isfield(opts, 'p0') && ~isempty(opts.p0)
  p = opts.p0(:); best = ref.window_cost(p, W);
else
  best = inf; p = [0; 0; 0];
  for th = th_g
    for dd = dd_g
      for d0 = d0_g
        c = ref.window_cost([th; d0; dd], W);
        if c < best, best = c; p = [th; d0; dd]; end
      end
    end
  end
end
info.grid_cost = best;

% --- 2. Gauss-Newton with a halving line search
cost = best;
h = [1e-4; 1e-4; 1e-6];                     % finite-difference steps
for it = 1:opt(opts, 'max_iter', 30)
  [~, ~, r] = ref.window_cost(p, W);
  J = zeros(numel(r), 3);
  for j = 1:3
    pj = p; pj(j) = pj(j) + h(j);
    [~, ~, rj] = ref.window_cost(pj, W);
    J(:, j) = (rj - r) / h(j);
  end
  dp = -(J' * J) \ (J' * r);
  step = 1; improved = false;
  for k = 1:20
    pt = p + step * dp;
    ct = ref.window_cost(pt, W);
    if ct < cost, improved = true; break; end
    step = step / 2;
  end
  if ~improved, break; end
  done = (cost - ct) < 1e-9 * cost;
  p = pt; cost = ct;
  if done, break; end
end
info.n_iter = it;

% --- 3. convention: ddelta >= 0, theta0 in [0, pi), delta0 in (-pi, pi]
if p(3) < 0
  p = [p(1) + pi/2; -p(2); -p(3)];
end
p(1) = mod(p(1), pi);
p(2) = angle(exp(1i * p(2)));
end

function v = opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
