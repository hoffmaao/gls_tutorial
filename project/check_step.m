function ok = check_step(k, impl)
%CHECK_STEP Test one step of your fabric estimator against the reference.
%
% check_step(k)           test YOUR function for step k (1-6)
% check_step(k, 'ref')    run the same test on the reference solution
%                         (instructors: confirms the test itself passes)
% ok = check_step(...)    also return true/false
%
%   1 coh_model         2 window_cost       3 fit_window
%   4 window_errors     5 fit_column        6 bootstrap_column
%
% Every test uses data whose answer is known, and tells you what is wrong
% when it fails - read the message before changing code.

if nargin < 2, impl = 'student'; end
here = fileparts(mfilename('fullpath'));
if isempty(which('apres.syntheticSite')) || isempty(which('ptt.fujitaModel'))  % set up the path only if
  addpath(fileparts(here));                        % it is not set up already,
  tutorial_setup();                                % so a student's own folder
end                                                % stays ahead of student/
if isempty(which('ref.coh_model'))
  addpath(fullfile(here, 'solutions'), '-end');    % the reference, as ref.*
end

f = handles(impl);
names = {'coh_model', 'window_cost', 'fit_window', 'window_errors', ...
         'fit_column', 'bootstrap_column'};
fprintf('\n--- step %d: %s (%s) ---\n', k, names{k}, impl);
try
  switch k
    case 1, ok = test_model(f);
    case 2, ok = test_cost(f);
    case 3, ok = test_fit(f);
    case 4, ok = test_errors(f);
    case 5, ok = test_column(f);
    case 6, ok = test_bootstrap(f);
    otherwise, error('check_step:k', 'steps are numbered 1 to 6');
  end
catch err
  if strcmp(err.identifier, 'project:todo')
    fprintf('NOT STARTED: %s\n', err.message);
  else
    fprintf('ERROR while running your code:\n  %s\n', err.message);
    if ~isempty(err.stack)
      fprintf('  at %s, line %d\n', err.stack(1).name, err.stack(1).line);
    end
  end
  ok = false;
end
if ok, fprintf('PASS\n'); else, fprintf('FAIL\n'); end
if nargout == 0, clear ok; end
end

% =========================================================================
function f = handles(impl)
if strcmp(impl, 'ref')
  f = struct('coh_model', @ref.coh_model, 'window_cost', @ref.window_cost, ...
    'fit_window', @ref.fit_window, 'window_errors', @ref.window_errors, ...
    'fit_column', @ref.fit_column, 'bootstrap_column', @ref.bootstrap_column);
else
  f = struct('coh_model', @coh_model, 'window_cost', @window_cost, ...
    'fit_window', @fit_window, 'window_errors', @window_errors, ...
    'fit_column', @fit_column, 'bootstrap_column', @bootstrap_column);
end
end

function [W, truth] = test_window(theta_deg, dlam, zc, seed)
% One 60 m window of a synthetic site with a uniform fabric below the firn:
% real coherence rows with their correlated covariance (apres.coherenceField),
% so a fit that ignores the whitening matrices gets the wrong answer.
Q = apres.calibratePhase(apres.syntheticSite(struct('top_m', 40, ...
  'theta', deg2rad(theta_deg), 'dlam', dlam), struct('seed', seed, 'z_max', 1500)));
F = apres.coherenceField(Q, struct('z_range', [zc - 40, zc + 40]));
r = abs(F.z - zc) <= 30;
W = struct('psi', F.psi, 'z', F.z_eff(r), 'zc', zc, 'C', F.C(r, :), ...
  'Wh', F.Wh(:, :, r), 'rank', F.rank(r));
truth = [deg2rad(theta_deg); NaN; dlam * 2*pi*300e6*0.034 / (sqrt(3.171)*0.299792458e9)];
end

function ok = test_model(f)
psi = deg2rad(0:10:170); z = (170:10:230)';
p = [deg2rad(35); 0.7; 0.02];
H = f.coh_model(psi, z, 200, p);
Hr = ref.coh_model(psi, z, 200, p);
ok = isequal(size(H), size(Hr));
if ~ok
  fprintf('H is %s but should be %s: one ROW per depth, one COLUMN per azimuth.\n', ...
    mat2str(size(H)), mat2str(size(Hr)));
  return
end
err = max(abs(H(:) - Hr(:)));
ok = err < 1e-10;
fprintf('largest difference from the reference: %.2g\n', err);
if ~ok
  if max(abs(H(:) - conj(Hr(:)))) < 1e-10
    fprintf('Your H is the complex CONJUGATE of the right one: check the sign of 1i*mu.*sin(delta).\n');
  elseif max(abs(abs(H(:)) - abs(Hr(:)))) < 1e-10
    fprintf('Magnitudes right, phases wrong: check mu = cos(2*(psi - theta0)).\n');
  else
    fprintf('Check delta = delta0 + ddelta*(z - zc), and the floor den = max(den, 0.05).\n');
  end
  return
end
% the magnitude of a single column is the same at every azimuth
fprintf('|H| at one depth ranges %.3f to %.3f over azimuth (should be one number).\n', ...
  min(abs(H(4,:))), max(abs(H(4,:))));
end

function ok = test_cost(f)
W = test_window(35, 0.08, 500, 1);
p = [deg2rad(40); 0.5; 0.01];
[c, g, r] = f.window_cost(p, W);
[cr, gr, rr] = ref.window_cost(p, W);
ok = true;
if ~isreal(g) || abs(g - gr) > 1e-8 * max(1, abs(gr))
  fprintf('g = %s, should be %.6f. Whiten data and model with Wh first, then sum over every entry.\n', num2str(g), gr);
  ok = false;
end
if ~isequal(size(r), size(rr))
  fprintf('r is %s but should be %s: one column, row by row, Wh_i * [Re; Im] of each row.\n', ...
    mat2str(size(r)), mat2str(size(rr)));
  ok = false;
elseif max(abs(r - rr)) > 1e-8 * max(1, max(abs(rr)))
  fprintf('r differs from the reference by up to %.2g. Whiten each row: Wh_i * [Re C_i, Im C_i]'' (pagemtimes(W.Wh, d)).\n', max(abs(r - rr)));
  ok = false;
end
if abs(c - cr) > 1e-8 * cr
  fprintf('cost = %.4f, should be %.4f (= r''*r).\n', c, cr);
  ok = false;
end
if ok
  fprintf('cost %.1f with %d independent numbers in the window.\n', c, sum(W.rank));
end
end

function ok = test_fit(f)
cases = [35 0.08 500 11; 172 0.12 400 12; 95 0.06 600 13];
ok = true;
for i = 1:size(cases, 1)
  [W, truth] = test_window(cases(i,1), cases(i,2), cases(i,3), cases(i,4));
  tic; [ph, c] = f.fit_window(W); t = toc;
  [pr, cr] = ref.fit_window(W);
  dth = rad2deg(abs(angle(exp(2i*(ph(1) - truth(1))))) / 2);
  fprintf('case %d: truth theta %5.1f deg ddelta %.4f | yours %5.1f deg %.4f | cost %.1f (ref %.1f) | %.1f s\n', ...
    i, cases(i,1), truth(3), rad2deg(ph(1)), ph(3), c, cr, t);
  if ph(3) < 0
    fprintf('  ddelta is negative: apply the convention (part C).\n'); ok = false;
  end
  if ph(1) < 0 || ph(1) >= pi
    fprintf('  theta0 must be wrapped into [0, pi).\n'); ok = false;
  end
  if c > cr * 1.001 + 1e-6
    fprintf('  your cost is higher than the reference: a false minimum, or Gauss-Newton stopped early.\n'); ok = false;
  end
  if dth > 1 || abs(ph(3) - pr(3)) > 0.02 * pr(3)
    fprintf('  too far from the reference answer (theta off the truth by %.2f deg).\n', dth); ok = false;
  end
end
end

function ok = test_errors(f)
W = test_window(35, 0.08, 500, 2);
ph = ref.fit_window(W);
E = f.window_errors(ph, W);
Er = ref.window_errors(ph, W);
ok = true;
if ~isfield(E, 'sigma') || numel(E.sigma) ~= 3
  fprintf('E.sigma should hold 3 numbers (theta0, delta0, ddelta).\n'); ok = false; return
end
rel = max(abs(E.sigma(:) - Er.sigma(:)) ./ Er.sigma(:));
fprintf('sigma: yours %s, reference %s\n', mat2str(E.sigma', 3), mat2str(Er.sigma', 3));
fprintf('chi2/dof: yours %.3f, reference %.3f (dof %d)\n', E.chi2_dof, Er.chi2_dof, Er.dof);
if abs(E.dof - Er.dof) > 0
  fprintf('  dof should be sum(W.rank) - 4: the independent numbers, not numel(r).\n'); ok = false;
end
if abs(E.chi2_dof - Er.chi2_dof) > 1e-6
  fprintf('  chi2_dof = cost / dof.\n'); ok = false;
end
if rel > 0.02
  fprintf('  sigma differs by %.0f%%: C_M = inv(J''*J), then sqrt(diag), then widen.\n', 100*rel); ok = false;
end
end

function ok = test_column(f)
[Q, truth] = apres.syntheticSite(struct('top_m', [0; 200], 'theta', deg2rad([30; 30]), ...
  'dlam', [0.05; 0.1]), struct('z_max', 450, 'seed', 4));
Q = apres.calibratePhase(Q);
F = apres.coherenceField(Q, struct('z_range', [20 420]));
out = f.fit_column(F);
outr = ref.fit_column(F);
ok = isfield(out, 'zw') && isequal(out.zw, outr.zw);
if ~ok, fprintf('Window centres differ from the reference: use zw = (z(1)+win/2 : step : z(end)-win/2).''\n'); return; end
dl_t = interp1(truth.top_m, truth.dlam, out.zw, 'previous', 'extrap');
disp(table(round(out.zw), round(rad2deg(out.theta0), 1), round(out.dlam, 3), round(dl_t, 3), out.ok, ...
  round(rad2deg(outr.theta0), 1), round(outr.dlam, 3), outr.ok, ...
  'VariableNames', {'z', 'theta', 'dlam', 'dlam_true', 'ok', 'ref_theta', 'ref_dlam', 'ref_ok'}))
if ~isequal(out.ok, outr.ok)
  fprintf('Your abstentions differ from the reference: check the two rules.\n'); ok = false;
end
good = out.ok & outr.ok;
if any(abs(out.dlam(good) - outr.dlam(good)) > 1e-3)
  fprintf('dlam differs from the reference: dlam = ddelta / grad_per_dlam.\n'); ok = false;
end
if ~isfield(out, 'theta0_fit') || any(isnan(out.theta0_fit))
  fprintf('theta0_fit is missing: keep the unmasked theta0 and dlam (step 6 needs them).\n'); ok = false;
end
if any(abs(out.sigma_dlam(good) - outr.sigma_dlam(good)) > 0.02 * outr.sigma_dlam(good))
  fprintf('sigma_dlam differs from the reference: sigma of ddelta divided by grad_per_dlam.\n'); ok = false;
end
end

function ok = test_bootstrap(f)
[Q, ~] = apres.syntheticSite(struct('top_m', [0; 200], 'theta', deg2rad([179.7; 179.7]), ...   % axis at the 0 = 180 wrap
  'dlam', [0.05; 0.1]), struct('z_max', 1500, 'seed', 4));
F = apres.coherenceField(apres.calibratePhase(Q), struct('z_range', [20 420]));
out = ref.fit_column(F);
o = struct('n_rep', 4);
tic; B = f.bootstrap_column(out, F, o); t = toc;
Br = ref.bootstrap_column(out, F, o);
k = out.ok;
disp(table(round(out.zw(k)), round(rad2deg(B.sd_theta0(k)), 4), round(rad2deg(Br.sd_theta0(k)), 4), ...
  round(B.sd_dlam(k), 5), round(Br.sd_dlam(k), 5), B.n_ok(k), ...
  'VariableNames', {'z', 'sd_theta', 'ref_sd_theta', 'sd_dlam', 'ref_sd_dlam', 'n_ok'}))
fprintf('%.0f s for %d windows x %d replicas\n', t, nnz(k), o.n_rep);
ok = isequal(B.n_ok, Br.n_ok);
if ~ok, fprintf('Replica counts differ: use seed seed0 + 1000*i + k and count only reported refits.\n'); return; end
rel = max(abs(B.sd_theta0(k) - Br.sd_theta0(k)) ./ Br.sd_theta0(k));
if ~(rel < 1e-6)
  ok = false;
  if any(B.sd_theta0(k) > 10 * Br.sd_theta0(k))
    fprintf('theta0 spread is far too large: fold each replica onto the window''s axis first.\n');
  else
    fprintf('theta0 spread differs by %.0f%%: std over the reported replicas only.\n', 100*rel);
  end
end
if any(abs(B.sd_dlam(k) - Br.sd_dlam(k)) > 1e-6 * Br.sd_dlam(k))
  fprintf('dlam spread differs from the reference.\n'); ok = false;
end
end
