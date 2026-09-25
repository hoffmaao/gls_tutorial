function out = fit_column(F, opts)
%FIT_COLUMN Step 6: the whole column.
%
% out = fit_column(F, opts)
%
% F: one site's coherence field (apres.coherenceField): z, psi, C, sigma,
% fc. For each window of win_m metres, every step_m metres:
%   1. W from the rows within win_m/2 of zw(i)
%   2. p = fit_window(W, opts); E = window_errors(p, W); Jk = window_jackknife(p, W)
%   3. dlam = ddelta / grad_per_dlam, with
%          grad_per_dlam = 2*pi*fc*deps / (n_ice*c)   (lesson 2's g)
%   4. abstain (theta0 and dlam = NaN) when
%        - ddelta * win_m < min_turn: too little phase change across the
%          window to locate the axis, and dlam is indistinguishable from 0
%        - chi2_dof > chi2_max: the model does not describe the data
%
% opts: .win_m (60), .step_m (60), .min_turn (0.5), .chi2_max (5)
%
% out, one entry per window centre zw [Nw x 1]: theta0, dlam, sigma_theta0,
%   sigma_dlam (step 4), jk_theta0, jk_dlam (step 5), delta0, ddelta,
%   gamma, chi2_dof, ok, grad_per_dlam
%
% Check:  check_step(6)

if nargin < 2, opts = struct(); end
win  = opt(opts, 'win_m', 60);
step = opt(opts, 'step_m', 60);
c = 0.299792458e9; deps = 0.034; n_ice = sqrt(3.171);
gpd = 2*pi*F.fc*deps / (n_ice*c);

zw = (F.z(1) + win/2 : step : F.z(end) - win/2)';
Nw = numel(zw);
out = struct('zw', zw, 'theta0', nan(Nw,1), 'dlam', nan(Nw,1), ...
  'sigma_theta0', nan(Nw,1), 'sigma_dlam', nan(Nw,1), ...
  'jk_theta0', nan(Nw,1), 'jk_dlam', nan(Nw,1), 'delta0', nan(Nw,1), ...
  'ddelta', nan(Nw,1), 'gamma', nan(Nw,1), 'chi2_dof', nan(Nw,1), ...
  'ok', false(Nw,1), 'grad_per_dlam', gpd);

for i = 1:Nw
  % TODO 6a: W, the fit, both error bars, dlam.
  % TODO 6b: out.ok(i) from the two rules.
  error('project:todo', 'fit_column is not written yet - see the TODOs in project/student/fit_column.m');
end
out.theta0(~out.ok) = NaN;
out.dlam(~out.ok) = NaN;
end

function v = opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
