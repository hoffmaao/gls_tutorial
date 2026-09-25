%% Project - Build a least-squares fabric estimator
%
% Six functions in project/student/, each with a header and TODOs, and a
% check:
%
%   step 1  coh_model         forward model                  lesson 6
%   step 2  window_cost       weighted misfit, closed-form g lessons 1-2
%   step 3  fit_window        grid search + Gauss-Newton     lesson 6
%   step 4  window_errors     linearized error bars, chi2    lessons 1, 4
%   step 5  window_jackknife  resampled error bars           lessons 3, 5
%   step 6  fit_column        the whole column               all
%
% Lessons 8 and 9 describe the data. Run this guide one section at a
% time; after editing a function, run check_step(k) in the Command Window
% until it prints PASS.

clear; close all;
if isempty(which('tutorial_setup'))
  here = fileparts(mfilename('fullpath'));
  if isempty(here), here = pwd; end
  addpath(fileparts(here));
end
tutorial_setup();
root = fullfile(fileparts(which('tutorial_setup')), 'data', 'GHOST24_Polarimetric_pRES_OZ');

%% Step 1 - forward model
% Edit project/student/coh_model.m, then check_step(1).
psi = deg2rad(0:2:178);
z = (0:5:300)';
H = coh_model(psi, z, 150, [deg2rad(30); 0; 0.03]);
figure('Name', 'Step 1', 'Position', [100 100 800 380]);
subplot(1,2,1); imagesc(rad2deg(psi), z, angle(H)); colorbar;
xlabel('\psi [deg]'); ylabel('depth [m]'); title('arg C');
subplot(1,2,2); imagesc(rad2deg(psi), z, abs(H)); colorbar; clim([0 1.05]);
xlabel('\psi [deg]'); title('|C|');
% Compare with the phase panel of lesson 9. |C| is 1 everywhere for a
% single uniform column.

%% Step 2 - misfit of a guess
% Edit project/student/window_cost.m, then check_step(2).
Q = apres.calibratePhase(apres.loadQuadpolSite(fullfile(root, 'PpRES_20240551_003')));
F = apres.coherenceField(Q, struct('z_range', [20 1200]));
rows = abs(F.z - 700) <= 30;
W = struct('psi', F.psi, 'z', F.z(rows), 'zc', 700, 'C', F.C(rows, :), 'sigma', F.sigma(rows, :));
for guess = [0 30 60 90 120 150]
  [c, g] = window_cost([deg2rad(guess); 0; 0.015], W);
  fprintf('theta0 = %3d deg: cost %8.0f, g = %.2f\n', guess, c, g);
end
% A negative g means the guessed pattern is the opposite of the data's.
% The best fit (step 3) has g > 0.

%% Step 3 - the best guess
% Edit project/student/fit_window.m, then check_step(3).
tic; [p, c] = fit_window(W); t = toc;
fprintf('best fit at 700 m: theta0 %.1f deg, ddelta %.4f rad/m, cost %.0f (%.1f s)\n', ...
  rad2deg(p(1)), p(3), c, t);

%% Step 4 - error bars and chi-square
% Edit project/student/window_errors.m, then check_step(4).
E = window_errors(p, W);
fprintf('theta0 = %.1f +- %.2f deg, chi2/dof = %.2f\n', rad2deg(p(1)), rad2deg(E.sigma(1)), E.chi2_dof);

%% Step 4, continued - test the error bars
% 20 synthetic sites with the same fabric and fresh speckle; fit the same
% window in each and compare the scatter with the claimed sigma. About a
% minute.
truth = struct('top_m', 0, 'theta', deg2rad(35), 'dlam', 0.08);
R = 20; th = zeros(R, 1); sg = zeros(R, 1); chi2 = zeros(R, 1);
dl = zeros(R, 1); sgd = zeros(R, 1);
for k = 1:R
  Qk = apres.calibratePhase(apres.syntheticSite(truth, struct('seed', 100 + k, 'z_max', 600)));
  Fk = apres.coherenceField(Qk, struct('z_range', [20 580]));
  rk = abs(Fk.z - 400) <= 30;
  Wk = struct('psi', Fk.psi, 'z', Fk.z(rk), 'zc', 400, 'C', Fk.C(rk, :), 'sigma', Fk.sigma(rk, :));
  pk = fit_window(Wk);
  Ek = window_errors(pk, Wk);
  th(k) = deg2rad(35) + angle(exp(2i*(pk(1) - deg2rad(35)))) / 2;   % fold: it is an axis
  sg(k) = Ek.sigma(1); chi2(k) = Ek.chi2_dof;
  dl(k) = pk(3); sgd(k) = Ek.sigma(3);
end
fprintf('theta0: measured scatter %.2f deg, claimed %.2f deg -> ratio %.2f\n', ...
  rad2deg(std(th)), rad2deg(median(sg)), median(sg) / std(th));
fprintf('ddelta: measured scatter %.5f, claimed %.5f -> ratio %.2f\n', ...
  std(dl), median(sgd), median(sgd) / std(dl));
fprintf('median chi2/dof %.2f\n', median(chi2));
% Both ratios are well below 1 while chi2/dof is below 1: the errors are
% correlated (lesson 3). All 18 azimuths come from the same 4 channels.

%% Step 5 - resampled error bars
% Edit project/student/window_jackknife.m, then check_step(5). Then repeat
% the test with the jackknife sigma:
sj = zeros(R, 1); sjd = zeros(R, 1);
for k = 1:R
  Qk = apres.calibratePhase(apres.syntheticSite(truth, struct('seed', 100 + k, 'z_max', 600)));
  Fk = apres.coherenceField(Qk, struct('z_range', [20 580]));
  rk = abs(Fk.z - 400) <= 30;
  Wk = struct('psi', Fk.psi, 'z', Fk.z(rk), 'zc', 400, 'C', Fk.C(rk, :), 'sigma', Fk.sigma(rk, :));
  Jk = window_jackknife(fit_window(Wk), Wk);
  sj(k) = Jk.sigma(1); sjd(k) = Jk.sigma(3);
end
fprintf('theta0: measured %.2f deg, jackknife %.2f deg -> ratio %.2f\n', ...
  rad2deg(std(th)), rad2deg(median(sj)), median(sj) / std(th));
fprintf('ddelta: measured %.5f, jackknife %.5f -> ratio %.2f\n', ...
  std(dl), median(sjd), median(sjd) / std(dl));
fprintf('theta0 mean error %.2f deg\n', rad2deg(mean(th) - deg2rad(35)));
% The jackknife makes the ddelta (so dlam) error bar honest. For theta0 it
% still falls short by 2-3x: speckle moves theta0 in a way shared across
% rows, which leaving one row out cannot see, and there is a depth-
% dependent bias of up to a few degrees. Quote dlam with its jackknife
% sigma and theta0 with a floor of about 3 degrees (open question 5).

%% Step 6 - the whole column, known answer
% Edit project/student/fit_column.m, then check_step(6).
% One axis (30 deg) throughout, contrast increasing layer by layer.
[Qs, tr] = apres.syntheticSite(struct('top_m', (0:250:1000)', ...
  'theta', deg2rad(30) * ones(5, 1), 'dlam', [0.05; 0.1; 0.15; 0.2; 0.25]), struct('seed', 7));
Fs = apres.coherenceField(apres.calibratePhase(Qs), struct('z_range', [20 1200]));
out = fit_column(Fs);
zz = [tr.top_m'; [tr.top_m(2:end)' 1200]];          % each layer's top and bottom
figure('Name', 'Step 6', 'Position', [100 100 800 420]);
subplot(1,2,1);
plot(repelem(rad2deg(tr.theta'), 2), zz(:), 'k', 'LineWidth', 1.5); hold on;
errorbar(rad2deg(out.theta0), out.zw, rad2deg(out.jk_theta0), 'horizontal', 'o');
set(gca, 'YDir', 'reverse'); xlim([10 50]); xlabel('\theta_0 [deg]'); ylabel('depth [m]'); grid on;
title('\theta_0(z)'); legend('truth', 'estimate', 'Location', 'southwest');
subplot(1,2,2);
plot(repelem(tr.dlam', 2), zz(:), 'k', 'LineWidth', 1.5); hold on;
errorbar(out.dlam, out.zw, out.jk_dlam, 'horizontal', 'o');
set(gca, 'YDir', 'reverse'); xlim([0 0.3]); xlabel('\Delta\lambda'); grid on; title('\Delta\lambda(z)');
% The top layer (dlam 0.05) turns the phase too little in 60 m, so those
% windows abstain. Windows across a layer boundary report a blend.
%
% Now rerun with the axis turning with depth:
%     'theta', deg2rad([20; 25; 30; 35; 40])
% dlam survives; theta0 fails below ~600 m. The model assumes one axis for
% the whole column above each window; when the axis turns, the coherence
% at depth carries every layer above it. Handling that needs a model of the
% full column, which is what ptt.fabricGLS is (lesson 7).

%% Four real sites
% The estimator on the GHOST sites beside ptt.ershadiFabric (Ershadi et
% al. 2022). theta0 is clockwise from the H antenna, which points north.
figure('Name', 'Finale', 'Position', [50 50 1400 500]);
for s = 1:4
  name = sprintf('PpRES_20240551_%03d', s);
  Q = apres.calibratePhase(apres.loadQuadpolSite(fullfile(root, name)));
  F = apres.coherenceField(Q, struct('z_range', [20 1600]));
  good = F.z < F.z(find(F.snr_db < 12, 1));         % stop where the signal fades
  F.z = F.z(good); F.C = F.C(good, :); F.sigma = F.sigma(good, :);
  out = fit_column(F);
  Er = ptt.ershadiFabric(struct('hh', Q.hh, 'vv', Q.vv, 'hv', Q.hv, 'vh', Q.vh), Q.z, ...
    struct('fc', Q.fc, 'win_m', 10, 'grad_win_m', 30, 'coh_min', 0.3));
  subplot(2, 4, s);
  plot(mod(rad2deg(Er.theta), 180), Q.z, '.', 'Color', [0.75 0.75 0.75], 'MarkerSize', 2); hold on;
  errorbar(rad2deg(out.theta0), out.zw, rad2deg(out.jk_theta0), 'horizontal', 'o', 'LineWidth', 1.2);
  set(gca, 'YDir', 'reverse'); xlim([0 180]); ylim([0 1200]); grid on;
  title(name, 'Interpreter', 'none'); xlabel('\theta_0 [deg from north]');
  if s == 1, ylabel('depth [m]'); end
  subplot(2, 4, 4 + s);
  plot(Er.dlam, Q.z, '.', 'Color', [0.75 0.75 0.75], 'MarkerSize', 2); hold on;
  errorbar(out.dlam, out.zw, out.jk_dlam, 'horizontal', 'o', 'LineWidth', 1.2);
  set(gca, 'YDir', 'reverse'); xlim([0 0.5]); ylim([0 1200]); grid on;
  xlabel('\Delta\lambda'); if s == 1, ylabel('depth [m]'); end
  fprintf('%s: %d of %d windows reported, median chi2/dof %.1f\n', name, nnz(out.ok), numel(out.ok), median(out.chi2_dof));
end
% Grey: ptt.ershadiFabric. Orange: the estimator, with jackknife error bars.

%% OPEN QUESTIONS
% 1. Where do the two methods disagree, and which is better supported?
%    Use chi2_dof and the jackknife sigmas.
% 2. Some windows abstain because chi2_dof > 5. Plot chi2_dof against
%    delta0. Where is the model's denominator smallest?
% 3. Sites 002 and 004 have strong fabric. Try win_m = 40 and 80.
% 4. apres.calibratePhase chose between two branches that mirror theta0
%    about north. Flip it (add pi to a and b) and rerun one site. What
%    field information would decide it?
% 5. The Monte Carlo measured how far the theta0 jackknife falls short.
%    Calibrate an inflation factor on synthetics matched to each site
%    (same n_looks, similar dlam) and apply it. What does that assume?
