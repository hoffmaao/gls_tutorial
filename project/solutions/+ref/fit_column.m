function out = fit_column(F, opts)
%FIT_COLUMN Reference solution, step 5: slide the window down the column.
%
% out = ref.fit_column(F, opts)
%
% F is what apres.coherenceField returns (z, z_eff, psi, C, Wh, rank, fc).
% Rows are assigned to windows by z and modelled at z_eff. For
% each window of win_m metres, stepped by step_m: fit (step 3), error bars
% (step 4), then convert the phase rate to the fabric contrast
%     dlam = ddelta / grad_per_dlam,  grad_per_dlam = 2 pi fc deps / (n c)
% (lesson 2's g, the constant ptt.quadpolFabricLS uses).
%
% Abstain (NaN) where the window cannot say anything:
%   - ddelta * win_m < min_turn: too little phase change to locate the
%     axis, and dlam is indistinguishable from zero;
%   - chi2_dof > chi2_max: the model does not describe the data (lesson 4).
% The model holds dlam and theta constant inside a window, so windows must
% be short where the fabric is strong; 60 m windows stepped by 60 m are
% also independent of each other.
%
% opts: .win_m (60), .step_m (60), .min_turn (0.5), .chi2_max (5), plus
%       anything ref.fit_window takes
%
% out fields per window centre zw [Nw x 1]: theta0 (rad, [0, pi)), dlam,
%   sigma_theta0, sigma_dlam, delta0, ddelta, gamma, chi2_dof, ok,
%   grad_per_dlam, win_m, and theta0_fit, dlam_fit (before abstention)

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
  'ok', false(Nw,1), 'grad_per_dlam', gpd, 'win_m', win);
for i = 1:Nw
  rows = abs(F.z - zw(i)) <= win/2;
  W = struct('psi', F.psi, 'z', F.z_eff(rows), 'zc', zw(i), ...
    'C', F.C(rows, :), 'Wh', F.Wh(:, :, rows), 'rank', F.rank(rows));
  p = ref.fit_window(W, opts);
  E = ref.window_errors(p, W);
  out.theta0(i) = p(1); out.delta0(i) = p(2); out.ddelta(i) = p(3);
  out.dlam(i) = p(3) / gpd;
  out.sigma_theta0(i) = E.sigma(1);
  out.sigma_dlam(i) = E.sigma(3) / gpd;
  out.gamma(i) = E.g;
  out.chi2_dof(i) = E.chi2_dof;
  out.ok(i) = p(3) * win >= opt(opts, 'min_turn', 0.5) && ...
              E.chi2_dof <= opt(opts, 'chi2_max', 5);
end
out.theta0_fit = out.theta0;          % unmasked, for step 6
out.dlam_fit = out.dlam;
out.theta0(~out.ok) = NaN;
out.dlam(~out.ok) = NaN;
end

function v = opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
