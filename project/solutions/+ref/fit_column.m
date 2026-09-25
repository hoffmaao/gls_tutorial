function out = fit_column(F, opts)
%FIT_COLUMN Reference solution, step 6: slide the window down the column.
%
% out = ref.fit_column(F, opts)
%
% F is what apres.observables returns (z, psi, C, sigma, fc). For each
% window of win_m metres, stepped by step_m: fit (step 3), error bars
% (step 4) and resampled error bars (step 5), then convert the phase
% rate to the fabric contrast
%     dlam = ddelta / grad_per_dlam,  grad_per_dlam = 2 pi fc deps / (n c)
% (lesson 2's g, the same constant ptt.quadpolFabricLS uses).
%
% ABSTAIN (report NaN) where the window cannot say anything:
%   - the phase barely rotates across it (ddelta * win_m < min_turn rad):
%     then there is no pattern whose position gives theta0 (lesson 6,
%     TRY THIS 1), and dlam is indistinguishable from zero;
% WINDOW LENGTH. The model holds dlam and theta constant inside a window,
% so the window must be short where the fabric is strong (at dlam 0.3 the
% phase turns a full cycle in ~170 m at 300 MHz). 60 m windows stepped by
% 60 m are also independent of each other.
%   - chi2_dof > chi2_max: the model does not describe the data (lesson 4).
%
% opts: .win_m (60), .step_m (60), .min_turn (0.5), .chi2_max (5), plus
%       anything ref.fit_window takes
%
% out fields per window centre zw [Nw x 1]: theta0 (rad, [0, pi)), dlam,
%   sigma_theta0, sigma_dlam (linearized, step 4), jk_theta0, jk_dlam
%   (jackknife, step 5 - the ones to quote), delta0, ddelta, gamma, chi2_dof,
%   ok (logical: not abstained)

if nargin < 2, opts = struct(); end
win  = opt(opts, 'win_m', 60);
step = opt(opts, 'step_m', 60);
c = 0.299792458e9; deps = 0.034; n_ice = sqrt(3.171);
gpd = 2*pi*F.fc*deps / (n_ice*c);                 % rad/m per unit dlam

zw = (F.z(1) + win/2 : step : F.z(end) - win/2)';
Nw = numel(zw);
out = struct('zw', zw, 'theta0', nan(Nw,1), 'dlam', nan(Nw,1), ...
  'sigma_theta0', nan(Nw,1), 'sigma_dlam', nan(Nw,1), 'delta0', nan(Nw,1), ...
  'ddelta', nan(Nw,1), 'gamma', nan(Nw,1), 'chi2_dof', nan(Nw,1), ...
  'ok', false(Nw,1), 'grad_per_dlam', gpd, ...
  'jk_theta0', nan(Nw,1), 'jk_dlam', nan(Nw,1));
for i = 1:Nw
  rows = abs(F.z - zw(i)) <= win/2;
  W = struct('psi', F.psi, 'z', F.z(rows), 'zc', zw(i), ...
    'C', F.C(rows, :), 'sigma', F.sigma(rows, :));
  p = ref.fit_window(W, opts);
  E = ref.window_errors(p, W);
  Jk = ref.window_jackknife(p, W);
  out.theta0(i) = p(1); out.delta0(i) = p(2); out.ddelta(i) = p(3);
  out.dlam(i) = p(3) / gpd;
  out.sigma_theta0(i) = E.sigma(1);
  out.sigma_dlam(i) = E.sigma(3) / gpd;
  out.gamma(i) = E.g;
  out.jk_theta0(i) = Jk.sigma(1);
  out.jk_dlam(i) = Jk.sigma(3) / gpd;
  out.chi2_dof(i) = E.chi2_dof;
  out.ok(i) = p(3) * win >= opt(opts, 'min_turn', 0.5) && ...
              E.chi2_dof <= opt(opts, 'chi2_max', 5);
end
out.theta0(~out.ok) = NaN;
out.dlam(~out.ok) = NaN;
end

function v = opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
