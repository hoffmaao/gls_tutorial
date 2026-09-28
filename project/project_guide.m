%% Project - Build a least-squares fabric estimator
%
% Six functions in project/student/, each with a header and TODOs, and a
% check:
%
%   step 1  coh_model         forward model                    lesson 6
%   step 2  window_cost       GLS misfit, closed-form g        lessons 1-3
%   step 3  fit_window        grid search + Gauss-Newton       lesson 6
%   step 4  window_errors     formula error bars, chi2         lessons 1, 4
%   step 5  fit_column        the whole column                 all
%   step 6  bootstrap_column  error bars by simulated repeats  lesson 1's habit
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
% One 60 m window of a real site. Each row's 18 azimuths come from the
% same four channels, so they carry only a few independent numbers
% (W.rank); W.Wh turns them into independent, unit-variance ones.
Q = apres.calibratePhase(apres.loadQuadpolSite(fullfile(root, 'PpRES_20240551_003')));
F = apres.coherenceField(Q, struct('z_range', [20 1200]));
rows = abs(F.z - 700) <= 30;
W = struct('psi', F.psi, 'z', F.z_eff(rows), 'zc', 700, 'C', F.C(rows, :), ...
  'Wh', F.Wh(:, :, rows), 'rank', F.rank(rows));
fprintf('%d rows x %d azimuths = %d numbers, of which %d are independent\n', ...
  numel(W.z), numel(W.psi), 2*numel(W.z)*numel(W.psi), sum(W.rank));
for guess = [0 30 60 90 120 150]
  [c, g] = window_cost([deg2rad(guess); 0; 0.015], W);
  fprintf('theta0 = %3d deg: cost %10.0f, g = %.2f\n', guess, c, g);
end
% g = 0 means the guessed pattern points the wrong way: no scale helps.

%% Step 3 - the best guess
% Edit project/student/fit_window.m, then check_step(3).
tic; [p, c] = fit_window(W); t = toc;
fprintf('best fit at 700 m: theta0 %.1f deg, ddelta %.4f rad/m, cost %.0f (%.1f s)\n', ...
  rad2deg(p(1)), p(3), c, t);

%% Step 4 - formula error bars and chi-square
% Edit project/student/window_errors.m, then check_step(4).
E = window_errors(p, W);
fprintf('theta0 = %.1f +- %.2f deg, chi2/dof = %.2f\n', rad2deg(p(1)), rad2deg(E.sigma(1)), E.chi2_dof);

%% Step 4, continued - test the formula
% 20 synthetic sites with the same fabric and fresh speckle; fit the same
% window in each and compare the scatter with the formula. About a minute.
truth = struct('top_m', 40, 'theta', deg2rad(35), 'dlam', 0.08);
R = 20; th = zeros(R, 1); sg = zeros(R, 1); dl = th; sgd = th; chi2 = th;
for k = 1:R
  Qk = apres.calibratePhase(apres.syntheticSite(truth, struct('seed', 100 + k, 'z_max', 1500)));
  Fk = apres.coherenceField(Qk, struct('z_range', [460 540]));
  rk = abs(Fk.z - 500) <= 30;
  Wk = struct('psi', Fk.psi, 'z', Fk.z_eff(rk), 'zc', 500, 'C', Fk.C(rk, :), ...
    'Wh', Fk.Wh(:, :, rk), 'rank', Fk.rank(rk));
  pk = fit_window(Wk);
  Ek = window_errors(pk, Wk);
  th(k) = deg2rad(35) + angle(exp(2i*(pk(1) - deg2rad(35)))) / 2;   % fold: it is an axis
  sg(k) = Ek.sigma(1); dl(k) = pk(3); sgd(k) = Ek.sigma(3); chi2(k) = Ek.chi2_dof;
end
fprintf('theta0: scatter %.3f deg, formula %.3f deg -> ratio %.1f\n', ...
  rad2deg(std(th)), rad2deg(median(sg)), median(sg) / std(th));
fprintf('ddelta: scatter %.1e, formula %.1e -> ratio %.1f\n', std(dl), median(sgd), median(sgd) / std(dl));
fprintf('theta0 mean error %.3f deg, median chi2/dof %.2f\n', rad2deg(mean(th)) - 35, median(chi2));
% A ratio above 1 means the formula is cautious. It is, deliberately:
% apres.coherenceField limits how precisely any one direction of the data
% can be trusted (cond_floor), because the exact covariance claims some
% directions to be known almost perfectly and small model errors there
% then dominate. Cautious is safe but not accurate; step 6 gets closer.

%% Step 5 - the whole column, known answer
% Edit project/student/fit_column.m, then check_step(5).
% One axis (30 deg) throughout, contrast increasing layer by layer.
[Qs, tr] = apres.syntheticSite(struct('top_m', (0:250:1000)', ...
  'theta', deg2rad(30) * ones(5, 1), 'dlam', [0.05; 0.1; 0.15; 0.2; 0.25]), ...
  struct('seed', 7, 'z_max', 3000));
Fs = apres.coherenceField(apres.calibratePhase(Qs), struct('z_range', [20 1200]));
out = fit_column(Fs);
zz = [tr.top_m'; [tr.top_m(2:end)' 1200]];          % each layer's top and bottom
figure('Name', 'Step 5', 'Position', [100 100 800 420]);
subplot(1,2,1);
plot(repelem(rad2deg(tr.theta'), 2), zz(:), 'k', 'LineWidth', 1.5); hold on;
errorbar(rad2deg(out.theta0), out.zw, rad2deg(out.sigma_theta0), 'horizontal', 'o');
set(gca, 'YDir', 'reverse'); xlim([25 35]); xlabel('\theta_0 [deg]'); ylabel('depth [m]'); grid on;
title('\theta_0(z)'); legend('truth', 'estimate', 'Location', 'southwest');
subplot(1,2,2);
plot(repelem(tr.dlam', 2), zz(:), 'k', 'LineWidth', 1.5); hold on;
errorbar(out.dlam, out.zw, out.sigma_dlam, 'horizontal', 'o');
set(gca, 'YDir', 'reverse'); xlim([0 0.3]); xlabel('\Delta\lambda'); grid on; title('\Delta\lambda(z)');
% Error bars: the formula (step 4). The top layer (dlam 0.05) turns the
% phase too little in 60 m, so those windows abstain; windows where the
% phase passes pi misfit and abstain too.
%
% Now rerun with the axis turning with depth:
%     'theta', deg2rad([20; 25; 30; 35; 40])
% Most windows abstain and the rest get theta0 wrong: the model assumes
% one axis for the whole column above each window. A turning axis needs a
% model of the full column, which is what ptt.fabricGLS is (lesson 7).

%% Step 6 - error bars by simulated repeats
% Edit project/student/bootstrap_column.m, then check_step(6).
% Repeat each window's measurement 20 times in simulation, and compare
% with the true scatter over 20 independent versions of the whole site
% (possible here only because the site is synthetic). About two minutes.
B = bootstrap_column(out, Fs, struct('n_rep', 20));
TH = nan(numel(out.zw), 20); DL = TH;
for k = 1:20
  Qk = apres.syntheticSite(tr, struct('seed', 500 + k, 'z_max', 3000));
  ok = fit_column(apres.coherenceField(apres.calibratePhase(Qk), struct('z_range', [20 1200])));
  TH(:, k) = deg2rad(30) + angle(exp(2i*(ok.theta0 - deg2rad(30)))) / 2;
  DL(:, k) = ok.dlam;
end
mc_th = std(TH, 0, 2, 'omitnan'); mc_dl = std(DL, 0, 2, 'omitnan');
k = out.ok;
disp(table(round(out.zw(k)), round(rad2deg(out.sigma_theta0(k)), 3), round(rad2deg(B.sd_theta0(k)), 3), ...
  round(rad2deg(mc_th(k)), 3), round(out.sigma_dlam(k), 4), round(B.sd_dlam(k), 4), round(mc_dl(k), 4), ...
  'VariableNames', {'z', 'theta_formula', 'theta_repeats', 'theta_true', 'dlam_formula', 'dlam_repeats', 'dlam_true'}))
figure('Name', 'Step 6', 'Position', [100 100 800 420]);
subplot(1,2,1);
semilogx(rad2deg(out.sigma_theta0(k)), out.zw(k), 's', rad2deg(B.sd_theta0(k)), out.zw(k), 'o', ...
  rad2deg(mc_th(k)), out.zw(k), 'k.', 'MarkerSize', 8, 'LineWidth', 1.2);
set(gca, 'YDir', 'reverse'); grid on; xlabel('\sigma of \theta_0 [deg]'); ylabel('depth [m]');
title('\theta_0 error bars'); legend('formula', 'simulated repeats', 'true scatter', 'Location', 'southeast');
subplot(1,2,2);
semilogx(out.sigma_dlam(k), out.zw(k), 's', B.sd_dlam(k), out.zw(k), 'o', mc_dl(k), out.zw(k), 'k.', ...
  'MarkerSize', 8, 'LineWidth', 1.2);
set(gca, 'YDir', 'reverse'); grid on; xlabel('\sigma of \Delta\lambda'); title('\Delta\lambda error bars');
% The repeats follow the true scatter window by window; the formula sits
% above both. On real data, only the formula and the repeats exist.

%% Four real sites
% The estimator on the GHOST sites beside ptt.ershadiFabric (Ershadi et
% al. 2022). Error bars: simulated repeats (10 per window). theta0 is
% clockwise from the H antenna, which points north. About a minute a site.
figure('Name', 'Finale', 'Position', [50 50 1400 500]);
for s = 1:4
  name = sprintf('PpRES_20240551_%03d', s);
  Q = apres.calibratePhase(apres.loadQuadpolSite(fullfile(root, name)));
  F = apres.coherenceField(Q, struct('z_range', [20 1600]));
  g = F.z < F.z(find(F.snr_db < 12, 1));             % stop where the signal fades
  for fld = {'z', 'z_eff', 'snr_db'}, F.(fld{1}) = F.(fld{1})(g); end
  F.C = F.C(g, :); F.Wh = F.Wh(:, :, g); F.rank = F.rank(g);
  out = fit_column(F);
  B = bootstrap_column(out, F, struct('n_rep', 10));
  Er = ptt.ershadiFabric(struct('hh', Q.hh, 'vv', Q.vv, 'hv', Q.hv, 'vh', Q.vh), Q.z, ...
    struct('fc', Q.fc, 'win_m', 10, 'grad_win_m', 30, 'coh_min', 0.3));
  subplot(2, 4, s);
  plot(mod(rad2deg(Er.theta), 180), Q.z, '.', 'Color', [0.75 0.75 0.75], 'MarkerSize', 2); hold on;
  errorbar(rad2deg(out.theta0), out.zw, rad2deg(B.sd_theta0), 'horizontal', 'o', 'LineWidth', 1.2);
  set(gca, 'YDir', 'reverse'); xlim([0 180]); ylim([0 1200]); grid on;
  title(name, 'Interpreter', 'none'); xlabel('\theta_0 [deg from north]');
  if s == 1, ylabel('depth [m]'); end
  subplot(2, 4, 4 + s);
  plot(Er.dlam, Q.z, '.', 'Color', [0.75 0.75 0.75], 'MarkerSize', 2); hold on;
  errorbar(out.dlam, out.zw, B.sd_dlam, 'horizontal', 'o', 'LineWidth', 1.2);
  set(gca, 'YDir', 'reverse'); xlim([0 0.5]); ylim([0 1200]); grid on;
  xlabel('\Delta\lambda'); if s == 1, ylabel('depth [m]'); end
  fprintf('%s: %d of %d windows reported, median chi2/dof %.1f\n', name, nnz(out.ok), numel(out.ok), median(out.chi2_dof));
end
% Grey: ptt.ershadiFabric. Orange: the estimator, error bars from repeats.

%% OPEN QUESTIONS
% 1. Where do the two methods disagree, and which is better supported?
%    Use chi2_dof and the error bars.
% 2. Some windows abstain because chi2_dof > 5. Plot chi2_dof against
%    delta0. Where is the model's denominator smallest?
% 3. Change cond_floor in apres.coherenceField (opts.cond_floor, default
%    0.05) to 0.01 and 0.2. How do the formula error bars, the step-4
%    ratios and the number of reported real windows change?
% 4. apres.calibratePhase chose between two branches that mirror theta0
%    about north. Flip it (add pi to a and b) and rerun one site. What
%    field information would decide it?
% 5. The repeats assume the fitted window is the truth, with the same
%    looks, signal and VV coherence as the data. Which of those could be
%    wrong on real ice, and in which direction would it bias the error bars?
