%% Lesson 6 - Nonlinear least squares: when the model isn't a matrix
%
% At one depth the four channels let us synthesize the radar at any
% antenna azimuth psi. The HH-VV coherence C(psi) then depends on
%
%   theta0 - azimuth of the fabric axis
%   delta  - phase difference the two polarizations have built up
%
% through (ptt.quadpolFabricLS), with mu = cos 2(psi - theta0):
%
%   num = (1-mu^2)/2 + ((1+mu^2)/2) cos(delta) + i mu sin(delta)
%   den = 1 - ((1-mu^2)/2) (1 - cos(delta))
%   C   = gamma * num / den          (gamma <= 1: coherence loss)
%
% No matrix G gives C = G*m. The tools instead: grid search
% (quadpolFabricLS), Gauss-Newton and multi-start (fabricGLS). The price:
% error bars become approximate.

clear; close all; rng(6);

%% 1. The model (coh_model, at the end of this file), p = [theta0; delta; gamma]
psi   = deg2rad(0:2:178)';           % 90 synthetic azimuths
p_true = [deg2rad(35); 1.2; 0.8];
C_true = coh_model(psi, p_true);

figure('Name', 'Lesson 6a', 'Position', [100 100 1000 340]);
subplot(1,2,1); hold on;
for th = [20 35 80]
  plot(rad2deg(psi), angle(coh_model(psi, [deg2rad(th); 1.2; 0.8])), 'LineWidth', 1.5);
end
grid on; xlabel('azimuth \psi [deg]'); ylabel('arg C [rad]');
xlim([0 180]); title('Model phase vs azimuth, \delta = 1.2'); legend('\theta_0 = 20', '35', '80');
subplot(1,2,2); hold on;
for dl = [0.4 1.2 2.4]
  plot(rad2deg(psi), angle(coh_model(psi, [deg2rad(35); dl; 0.8])), 'LineWidth', 1.5);
end
grid on; xlabel('azimuth \psi [deg]'); ylabel('arg C [rad]');
xlim([0 180]); title('Model phase vs azimuth, \theta_0 = 35'); legend('\delta = 0.4', '1.2', '2.4');

% The phase peaks at +delta on the fabric axis and bottoms at -delta 90 deg
% away. For a uniform column |C| = gamma at every azimuth:
[min(abs(C_true)), max(abs(C_true))]
% So the phase pattern carries the fabric. (|C| varies in real, layered
% columns; that is fabricGLS's business.)

%% 2. Noisy data, and a residual that never wraps
% Compare model and data as complex numbers, real and imaginary parts
% stacked. Phases wrap: 3.1 and -3.1 rad are nearly equal, yet
3.1 - (-3.1)
angle(exp(1i*(3.1 - (-3.1))))        % the wrapped difference: correct
% fabricGLS fits phase and wraps every residual (forgetting it gave
% chi2/dof ~ 2200). Fitting the complex field avoids the issue.
sigma = 0.08;                        % noise per real and imaginary part
d = C_true + sigma*(randn(size(psi)) + 1i*randn(size(psi)));

resid  = @(p) [real(d - coh_model(psi, p)); imag(d - coh_model(psi, p))];
misfit = @(p) sum(resid(p).^2) / sigma^2;     % chi-square (lesson 4)

%% 3. Tool 1: grid search, with gamma solved exactly at every node
% gamma multiplies the model, so for fixed theta0 and delta it is a
% one-unknown lesson-1 problem: gamma = sum(Re(conj(H).*d)) / sum(|H|^2),
% H being the model with gamma = 1. quadpolFabricLS does this too, so its
% grid covers only the truly nonlinear unknowns.
th_grid = deg2rad(0:1:179);
dl_grid = linspace(0, 2*pi, 181);
MIS = zeros(numel(dl_grid), numel(th_grid));
for a = 1:numel(th_grid)
  for b = 1:numel(dl_grid)
    H = coh_model(psi, [th_grid(a); dl_grid(b); 1]);
    gam = sum(real(conj(H).*d)) / sum(abs(H).^2);
    gam = min(max(gam, 0), 1.05);   % a coherence scale is never negative or >1
    MIS(b, a) = sum(abs(d - gam*H).^2) / sigma^2;
  end
end
[~, ib] = min(MIS(:));
[b0, a0] = ind2sub(size(MIS), ib);
fprintf('grid best: theta0 = %.0f deg, delta = %.2f   (truth %.0f deg, %.2f)\n', ...
  rad2deg(th_grid(a0)), dl_grid(b0), rad2deg(p_true(1)), p_true(2));
% If that differs from the truth, see the figure.

figure('Name', 'Lesson 6b', 'Position', [100 100 620 460]);
imagesc(rad2deg(th_grid), dl_grid, log10(MIS)); axis xy; colorbar; hold on;
plot(rad2deg(p_true(1)), p_true(2), 'wp', 'MarkerSize', 14, 'MarkerFaceColor', 'w');
plot(rad2deg(p_true(1)) + 90, 2*pi - p_true(2), 'wo', 'MarkerSize', 12, 'LineWidth', 1.5);
xlabel('\theta_0 [deg]'); ylabel('\delta [rad]');
title('log_{10} misfit over \theta_0 and \delta');
% (theta0, delta) and (theta0 + 90, 2*pi - delta) give the same field:
% rotating 90 deg swaps the axes and flips the phase's sign. No data
% separates them; a convention picks one (quadpolFabricLS: the axis along
% which the phase grows with depth). theta0 is known modulo 180 deg.

%% 4. Tool 2: Gauss-Newton
% Near a guess p, a small change dp moves the prediction by about J*dp,
% where the JACOBIAN J has one column per parameter: how every prediction
% moves when that parameter is nudged. Built by finite differences, as in
% fabricGLS's H_jac:
p_guess = [deg2rad(50); 1.0; 0.7];
J = jacobian_fd(resid, p_guess);
size(J)                   % 180 rows (90 real + 90 imag) x 3 parameters

% "residual = J*dp" is lesson 1 for the step:
%     dp = (J' J) \ (J' r)
% Take the step, recompute J, repeat. If a step raises the misfit, halve it
% (a line search; fabricGLS's H_solve halves up to 20 times).
[p_gn, hist_good] = gauss_newton(resid, misfit, p_guess);
fprintf('\nGN from a nearby start: theta0 = %.1f deg, delta = %.3f, gamma = %.3f in %d steps\n', ...
  rad2deg(p_gn(1)), p_gn(2), p_gn(3), numel(hist_good)-1);

%% 5. How it fails: a false bottom
p_bad = [deg2rad(120); 3.0; 0.7];
[p_gn_bad, hist_bad] = gauss_newton(resid, misfit, p_bad);
fprintf('GN from a far start:    theta0 = %.1f deg, delta = %.3f, misfit %.0f (vs %.0f)\n', ...
  rad2deg(mod(p_gn_bad(1), pi)), p_gn_bad(2), hist_bad(end), hist_good(end));
% Gauss-Newton only walks downhill. fabricGLS, started from its prior mean
% on perfect synthetic data, once returned a reversed dlam profile.

figure('Name', 'Lesson 6c', 'Position', [100 100 560 380]);
semilogy(0:numel(hist_good)-1, hist_good, 'o-', 'LineWidth', 1.5); hold on;
semilogy(0:numel(hist_bad)-1, hist_bad, 's-', 'LineWidth', 1.5);
yline(numel(d)*2 - 3, '--', 'expected at the truth: N - P');
xlabel('Gauss-Newton step'); ylabel('misfit \chi^2'); grid on;
legend('near start', 'far start', 'Location', 'northeast');
title('Misfit per Gauss-Newton step');

%% 6. Tool 3: multi-start
% Score a coarse grid of starts without iterating, then run Gauss-Newton
% from the best few and keep the lowest (fabricGLS's n_start).
th_try = deg2rad(0:30:150); dl_try = [0.5 1.5 2.5 3.5 4.5 5.5];
starts = []; score = [];
for th = th_try
  for dl = dl_try
    starts(:, end+1) = [th; dl; 0.8];                  %#ok<SAGROW>
    score(end+1) = misfit(starts(:, end));             %#ok<SAGROW>
  end
end
[~, order] = sort(score);
best_mis = inf;
for q = 1:3
  [pq, hq] = gauss_newton(resid, misfit, starts(:, order(q)));
  if hq(end) < best_mis, best_mis = hq(end); p_hat = pq; end
end
p_hat(1) = mod(p_hat(1), pi);
fprintf('\nmulti-start: theta0 = %.1f deg, delta = %.3f, gamma = %.3f\n', ...
  rad2deg(p_hat(1)), p_hat(2), p_hat(3));

%% 7. Linearized error bars
% At the answer, use J as if it were G:  C_M = inv(J' inv(C_d) J). Exact
% only if the model is linear over the range the noise moves the answer.
% Check by Monte Carlo at two noise levels.
for sigma_mc = [0.08 0.30]
  Jh = jacobian_fd(@(p) [real(coh_model(psi, p)); imag(coh_model(psi, p))], p_hat);
  C_M = inv(Jh' * Jh / sigma_mc^2);
  claimed = sqrt(diag(C_M));
  R = 300; est = nan(3, R);
  for k = 1:R
    dk = C_true + sigma_mc*(randn(size(psi)) + 1i*randn(size(psi)));
    rk = @(p) [real(dk - coh_model(psi, p)); imag(dk - coh_model(psi, p))];
    mk = @(p) sum(rk(p).^2);
    pk = gauss_newton(rk, mk, p_hat);            % start at the answer: best case
    pk(1) = p_true(1) + angle(exp(2i*(pk(1) - p_true(1))))/2;   % fold: it is an axis
    est(:, k) = pk;
  end
  measured = std(est, 0, 2);
  fprintf('\nnoise %.2f      claimed    measured   ratio\n', sigma_mc);
  fprintf('theta0 [deg] %8.2f   %8.2f   %.2f\n', rad2deg(claimed(1)), rad2deg(measured(1)), claimed(1)/measured(1));
  fprintf('delta  [rad] %8.3f   %8.3f   %.2f\n', claimed(2), measured(2), claimed(2)/measured(2));
end
% Here the linearized error bars hold: three parameters, one well-sampled
% depth, fits started at the answer. On fabricGLS's 48-parameter layered
% problem its own test measures ~2x more scatter than claimed, so it calls
% its sigmas LOWER BOUNDS. Lesson 7 repeats that check.

%% TAKEAWAYS
% - Nonlinear least squares = lesson 1 repeated: linearize, step, repeat,
%   with a line search.
% - Periodic models have false bottoms and exact twins: grid or multi-start,
%   and know which ambiguities no data can break.
% - Linear nuisances (gamma) are solved in closed form inside the search.
% - Linearized error bars are claims; check them by Monte Carlo.
%
%% TRY THIS
% 0. In section 7, start each fit from the multi-start of section 6 and
%    raise the noise to 0.6. Do some land in a twin or false bottom?
% 1. Set p_true(2) = 0.1. Can theta0 still be found? What does
%    quadpolFabricLS do here? (See dlam_min_theta.)
% 2. Drop the line search. Does the near start still converge?
% 3. Add a pedestal: d = d + 0.1*sin(2*psi). What happens to theta0? This
%    is the instrument artefact quadpolFabricLS models.

%% ------------------------------------------------------------------------
function C = coh_model(psi, p)
% Single-column HH-VV coherence (ptt.quadpolFabricLS).
th = p(1); dl = p(2); gam = p(3);
mu  = cos(2*(psi - th));
num = (1 - mu.^2)/2 + ((1 + mu.^2)/2) .* cos(dl) + 1i * mu .* sin(dl);
den = 1 - ((1 - mu.^2)/2) .* (1 - cos(dl));
den = max(den, 0.05);      % 0/0 at delta = pi, mu = 0; the floor avoids NaN
C   = gam * num ./ den;
end

function J = jacobian_fd(f, p)
% Finite-difference Jacobian: nudge each parameter, see what moves.
f0 = f(p);
J = zeros(numel(f0), numel(p));
for j = 1:numel(p)
  h = 1e-6 * max(1, abs(p(j)));
  pj = p; pj(j) = pj(j) + h;
  J(:, j) = (f(pj) - f0) / h;
end
end

function [p, hist] = gauss_newton(resid, misfit, p)
% Gauss-Newton with a halving line search (fabricGLS's H_solve).
hist = misfit(p);
for it = 1:50
  r  = resid(p);
  J  = -jacobian_fd(resid, p);         % Jacobian of the model = -d(resid)/dp
  dp = (J' * J) \ (J' * r);
  step = 1; improved = false;
  for h = 1:20
    pt = p + step * dp;
    if misfit(pt) < hist(end), improved = true; break; end
    step = step / 2;
  end
  if ~improved, break; end
  p = pt;
  hist(end+1) = misfit(p);             %#ok<AGROW>
  if abs(hist(end-1) - hist(end)) < 1e-8 * hist(end), break; end
end
end
