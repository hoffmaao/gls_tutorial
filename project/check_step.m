function ok = check_step(k, impl)
%CHECK_STEP Test one step of your fabric estimator against the reference.
%
% check_step(k)           test YOUR function for step k (1-6)
% check_step(k, 'ref')    run the same test on the reference solution
%                         (instructors: confirms the test itself passes)
% ok = check_step(...)    also return true/false
%
%   1 coh_model         2 window_cost       3 fit_window
%   4 window_errors     5 window_jackknife  6 fit_column
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
         'window_jackknife', 'fit_column'};
fprintf('\n--- step %d: %s (%s) ---\n', k, names{k}, impl);
try
  switch k
    case 1, ok = test_model(f);
    case 2, ok = test_cost(f);
    case 3, ok = test_fit(f);
    case 4, ok = test_errors(f);
    case 5, ok = test_jackknife(f);
    case 6, ok = test_column(f);
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
    'window_jackknife', @ref.window_jackknife, 'fit_column', @ref.fit_column);
else
  f = struct('coh_model', @coh_model, 'window_cost', @window_cost, ...
    'fit_window', @fit_window, 'window_errors', @window_errors, ...
    'window_jackknife', @window_jackknife, 'fit_column', @fit_column);
end
end

function W = test_window(p, noise, seed)
% A window of data made from the model itself: 7 rows 10 m apart, 18
% azimuths, coherence 0.8, plus noise of the given size.
rng(seed);
psi = deg2rad(0:10:170);
z = (170:10:230)';
W = struct('psi', psi, 'z', z, 'zc', 200);
C0 = 0.8 * ref.coh_model(psi, z, 200, p);
W.sigma = noise * ones(size(C0));
W.C = C0 + noise * complex(randn(size(C0)), randn(size(C0)));
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
p = [deg2rad(35); 0.7; 0.02];
W = test_window(p, 0.05, 1);
[c, g, r] = f.window_cost(p, W);
[cr, gr, rr] = ref.window_cost(p, W);
ok = true;
if ~isreal(g) || abs(g - gr) > 1e-10
  fprintf('g = %s, should be %.6f. Use real(conj(H).*C) and sum over ALL entries.\n', num2str(g), gr);
  ok = false;
end
if ~isequal(size(r), size(rr))
  fprintf('r is %s but should be %s: one real COLUMN, real parts then imaginary parts.\n', ...
    mat2str(size(r)), mat2str(size(rr)));
  ok = false;
elseif max(abs(r - rr)) > 1e-10
  fprintf('r differs from the reference by up to %.2g. Divide by sigma; stack [real; imag].\n', max(abs(r - rr)));
  ok = false;
end
if abs(c - cr) > 1e-8 * cr
  fprintf('cost = %.4f, should be %.4f (= r''*r).\n', c, cr);
  ok = false;
end
if ok
  fprintf('cost %.1f for %d real residuals: about 1 per residual, as it should be at the truth.\n', c, numel(r));
end
end

function ok = test_fit(f)
cases = [35 0.7 0.02; 172 -2.5 0.035; 95 1.9 0.012];
ok = true;
for i = 1:size(cases, 1)
  p = [deg2rad(cases(i,1)); cases(i,2); cases(i,3)];
  W = test_window(p, 0.03, 10 + i);
  tic; [ph, c] = f.fit_window(W); t = toc;
  [~, cr] = ref.fit_window(W);
  dth = rad2deg(abs(angle(exp(2i*(ph(1) - p(1))))) / 2);
  fprintf('case %d: truth theta %5.1f deg ddelta %.3f | yours %5.1f deg %.4f | cost %.1f (ref %.1f) | %.1f s\n', ...
    i, cases(i,1), cases(i,3), rad2deg(ph(1)), ph(3), c, cr, t);
  if ph(3) < 0
    fprintf('  ddelta is negative: apply the convention (part C).\n'); ok = false;
  end
  if ph(1) < 0 || ph(1) >= pi
    fprintf('  theta0 must be wrapped into [0, pi).\n'); ok = false;
  end
  if c > cr * 1.001 + 1e-6
    fprintf('  your cost is higher than the reference: a false bottom, or Gauss-Newton stopped early.\n'); ok = false;
  end
  if dth > 1 || abs(ph(3) - p(3)) > 5e-4
    fprintf('  too far from the truth (theta off by %.2f deg).\n', dth); ok = false;
  end
end
end

function ok = test_errors(f)
p = [deg2rad(35); 0.7; 0.02];
W = test_window(p, 0.05, 2);
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
  fprintf('  dof should be numel(r) - 4.\n'); ok = false;
end
if abs(E.chi2_dof - Er.chi2_dof) > 1e-6
  fprintf('  chi2_dof = cost / dof.\n'); ok = false;
end
if rel > 0.02
  fprintf('  sigma differs by %.0f%%: C_M = inv(J''*J), then sqrt(diag), then widen.\n', 100*rel); ok = false;
end
end

function ok = test_jackknife(f)
% THE BRANCH TRAP. Put the ESTIMATED axis exactly on 0 = 180 deg, so the
% leave-one-out replicates land on both sides of the wrap: make the data,
% see where the estimate lands, then remake them (same noise) with the
% truth shifted by that much.
p = [0; -2.5; 0.035];
for it = 1:6                                      % walk the estimate onto the wrap
  W = test_window(p, 0.05, 3);
  ph = ref.fit_window(W);
  off = angle(exp(2i*ph(1))) / 2;                 % estimate's distance from 0 = 180
  if abs(off) < deg2rad(1e-3), break; end
  p(1) = p(1) - off;
end
J = f.window_jackknife(ph, W);
Jr = ref.window_jackknife(ph, W);
ok = isfield(J, 'sigma') && numel(J.sigma) == 3;
if ~ok, fprintf('J.sigma should hold 3 numbers.\n'); return; end
fprintf('jackknife sigma: yours %s, reference %s\n', mat2str(J.sigma', 3), mat2str(Jr.sigma', 3));
rel = max(abs(J.sigma(:) - Jr.sigma(:)) ./ Jr.sigma(:));
if rel > 0.05
  ok = false;
  if J.sigma(1) > 10 * Jr.sigma(1)
    fprintf('  theta0 sigma is huge: the true axis is at 0 = 180 deg, so replicates land on\n');
    fprintf('  near 0 deg. Fold each replicate onto the full answer''s branch first.\n');
  else
    fprintf('  differs by %.0f%%: drop ONE row per replicate, restart from p, scale by (n-1)/n.\n', 100*rel);
  end
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
if any(isnan(out.jk_dlam(outr.ok)))
  fprintf('jk_dlam is missing: store the jackknife sigma of ddelta divided by grad_per_dlam.\n'); ok = false;
end
end
