function B = bootstrap_column(out, F, opts)
%BOOTSTRAP_COLUMN Step 6: check the error bars by repeating each window's
%measurement in simulation (a parametric bootstrap).
%
% B = bootstrap_column(out, F, opts)
%
% out is your fit_column's result for the coherence field F. For every
% reported window: if its fitted values were the truth, how much would the
% answer change from one survey to the next? Build synthetic sites that
% reproduce that window, refit the window in each, and measure the spread.
% Each replica has the window's theta0 from the firn down, the window's
% dlam around the window, and an upper layer that puts the window at its
% fitted phase delta0 (the same point in the phase cycle). Windows centred
% shallower than about 115 m leave no room for that upper layer, so they
% are skipped (NaN spread, n_ok 0). Each replica is matched to the data:
%   signal   decay_m from a straight-line fit of F.snr_db (dB) against
%            depth, anchored to the window's own SNR
%   VV       decorrelation so the replica's fitted g matches the data's
%            (a pilot replica measures the g the replica gives without it)
%   looks    speckle_width chosen so the synthetic n_looks matches F's
% Replicate k of window i uses seed seed0 + 1000*i + k; the pilot uses
% seed0 + 1000*i (reproducible).
%
% The replica construction is written for you below (read it: every
% choice is there to make the replica match the real window). Your part is
% the experiment: repeat, refit, measure the spread.
%
% opts: .n_rep (20), .seed0 (1000), plus anything fit_column takes
%
% Check:  check_step(6)
%
% B fields per window: sd_theta0, sd_dlam (spread of the refits; NaN for
%   abstained or skipped shallow windows or fewer than 3 reported refits),
%   n_ok (refits that reported), theta0_reps, dlam_reps [Nw x n_rep], site_opts

if nargin < 3, opts = struct(); end
R   = opt(opts, 'n_rep', 20);
s0  = opt(opts, 'seed0', 1000);
win = out.win_m;
gpd = out.grad_per_dlam;
firn = 40;
Nw  = numel(out.zw);

% --- matched signal and looks, shared by all windows
cf = polyfit(F.z(:), F.snr_db(:), 1);
base = struct('snr_surface_db', cf(2), ...
  'decay_m', min(max(-(10/log(10)) / min(cf(1), -eps), 20), 5000), ...
  'speckle_width', matchWidth(F.n_looks, F.win_m, F.dz), ...
  'dz', F.dz, 'fc', F.fc, 'firn_m', firn);

B = struct('sd_theta0', nan(Nw,1), 'sd_dlam', nan(Nw,1), 'n_ok', zeros(Nw,1), ...
  'theta0_reps', nan(Nw,R), 'dlam_reps', nan(Nw,R), 'site_opts', base);
fold = @(t, t0) t0 + angle(exp(2i*(t - t0))) / 2;     % theta0 is an axis

for i = find(out.ok(:)).'
  zc = out.zw(i); th = out.theta0(i); dd = out.ddelta(i);
  % Replica column: the window's axis from the firn down, so there is
  % cross-polarized signal for apres.calibratePhase (it needs HV and VH
  % between 40 and 800 m). The window sits in a lower layer with the
  % window's dlam; the upper layer's dlam is chosen so the phase at the
  % window equals the fitted delta0. Too shallow for that layer: skip.
  top = zc - win/2 - 25;
  if top <= firn + 20, continue; end
  need = mod(out.delta0(i) - dd * (zc - top), 2*pi) + 2*pi*(0:50);
  dl_up = need / (gpd * (top - firn));
  dl_up = dl_up(find(dl_up >= 0.02, 1));
  truth = struct('top_m', [firn; top], 'theta', [th; th], 'dlam', [dl_up; dd / gpd]);
  so = base;
  so.z_max = max(zc + win/2 + 40, 1500);   % a full-length record for coregistration
  % signal: the fitted decay, anchored to this window's measured SNR
  snr_db = interp1(F.z, F.snr_db, zc, 'linear', 'extrap');
  so.snr_surface_db = snr_db + (10/log(10)) * zc / so.decay_m;
  % VV decorrelation: the fitted g is below 1 partly for reasons the
  % replica reproduces by itself (receiver noise, the phase turning within
  % each row). A pilot replica without decorrelation measures that part;
  % impose only the remainder, g_data / g_pilot.
  so.vv_coherence = 1;
  so.seed = s0 + 1000*i;
  g_pilot = refitWindow(truth, so, zc, win, F.win_m, opts);
  if isfinite(g_pilot.g) && g_pilot.g > 0
    so.vv_coherence = min(max(out.gamma(i) / g_pilot.g, 0.05), 1);
  end
  % TODO 6a: for k = 1:R, set so.seed = s0 + 1000*i + k, call
  %          rk = refitWindow(truth, so, zc, win, F.win_m, opts), and when
  %          rk.ok store fold(rk.theta0, th) in B.theta0_reps(i, k) and
  %          rk.dlam in B.dlam_reps(i, k). (fold: theta0 is an axis.)
  % TODO 6b: B.n_ok(i) = number of finite replicas; if at least 3,
  %          B.sd_theta0(i) and B.sd_dlam(i) = their standard deviations
  %          (std(..., 'omitnan')).
  error('project:todo', 'bootstrap_column is not written yet - see TODO 6a-6b');
end
end

function r = refitWindow(truth, so, zc, win, cwin, opts)
% One replica: make the site, build its coherence field around the window,
% and refit the window nearest zc.
Q = apres.calibratePhase(apres.syntheticSite(truth, so));
Fk = apres.coherenceField(Q, struct('z_range', [zc - win/2 - 5, zc + win/2 + 15], 'win_m', cwin));
opts.win_m = win; opts.step_m = win;
o1 = fit_column(Fk, opts);
[~, j] = min(abs(o1.zw - zc));
r = struct('ok', ~isempty(j) && o1.ok(j), 'theta0', o1.theta0_fit(j), ...
  'dlam', o1.dlam_fit(j), 'g', o1.gamma(j));
if isempty(j), r.g = NaN; end
end

function w = matchWidth(n_target, win_m, dz)
% Gaussian speckle width [bins] whose window of win_m metres holds n_target
% independent looks (the formula apres.observables uses).
nr = 2*floor(win_m / dz / 2) + 1;
lo = 0.3; hi = 20;
for it = 1:40
  w = (lo + hi) / 2;
  if nLooks(apres.syntheticSite_kernel(w), nr) > n_target, lo = w; else, hi = w; end
end
end

function n = nLooks(k, nr)
c = conv(k, flipud(k)); c = c / max(c);
m = (numel(c) + 1) / 2;
rho = zeros(nr, 1);
L = min(nr - 1, m - 1);
rho(1:L+1) = c(m:m+L);
w = 1 - (0:nr-1)' / nr;
n = nr / (1 + 2 * sum(w(2:end) .* rho(2:end).^2));
end

function v = opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
