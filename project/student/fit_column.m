function out = fit_column(F, opts)
%FIT_COLUMN Step 5: the whole column.
%
% out = fit_column(F, opts)
%
% F: one site's coherence field (apres.coherenceField). For each window of
% win_m metres, every step_m metres:
%   1. rows = those whose centre F.z lies within win_m/2 of zw(i);
%      W = struct('psi', F.psi, 'z', F.z_eff(rows), 'zc', zw(i),
%                 'C', F.C(rows,:), 'Wh', F.Wh(:,:,rows), 'rank', F.rank(rows))
%      z_eff, not z: each row's coherence is a brightness-weighted average
%      over its 10 m, and belongs to that average depth.
%   2. p = fit_window(W, opts); E = window_errors(p, W)
%   3. dlam = ddelta / grad_per_dlam, with
%          grad_per_dlam = 2*pi*fc*deps / (n_ice*c)   (lesson 2's g)
%   4. abstain (theta0 and dlam = NaN) when
%        - ddelta * win_m < min_turn: too little phase change across the
%          window to locate the axis, and dlam is indistinguishable from 0
%        - chi2_dof > chi2_max: the model does not describe the data
%      Keep the unmasked values in theta0_fit and dlam_fit (step 6 uses them).
%
% opts: .win_m (60), .step_m (60), .min_turn (0.5), .chi2_max (5)
%
% out, one entry per window centre zw [Nw x 1]: theta0, dlam, theta0_fit,
%   dlam_fit, sigma_theta0, sigma_dlam, delta0, ddelta, gamma, chi2_dof,
%   ok, grad_per_dlam, win_m
%
% Check:  check_step(5)

if nargin < 2, opts = struct(); end
win  = opt(opts, 'win_m', 60);
step = opt(opts, 'step_m', 60);
c = 0.299792458e9; deps = 0.034; n_ice = sqrt(3.171);
gpd = 2*pi*F.fc*deps / (n_ice*c);

zw = (F.z(1) + win/2 : step : F.z(end) - win/2)';
Nw = numel(zw);
out = struct('zw', zw, 'theta0', nan(Nw,1), 'dlam', nan(Nw,1), ...
  'sigma_theta0', nan(Nw,1), 'sigma_dlam', nan(Nw,1), 'delta0', nan(Nw,1), ...
  'ddelta', nan(Nw,1), 'gamma', nan(Nw,1), 'chi2_dof', nan(Nw,1), ...
  'ok', false(Nw,1), 'grad_per_dlam', gpd, 'win_m', win);

for i = 1:Nw
  % TODO 5a: W, the fit, the error bars, dlam; store every field.
  % TODO 5b: out.ok(i) from the two rules.
  error('project:todo', 'fit_column is not written yet - see the TODOs in project/student/fit_column.m');
end
out.theta0_fit = out.theta0;
out.dlam_fit = out.dlam;
out.theta0(~out.ok) = NaN;
out.dlam(~out.ok) = NaN;
end

function v = opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
