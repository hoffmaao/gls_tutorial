%% Lesson 7 - Reading ptt.fabricGLS
%
% Keep vendor/+ptt/fabricGLS.m open beside this file. Each block below
% names the code it refers to and the lesson it came from. Then the
% function runs on a synthetic column and its error bars are checked as
% its own test does.
%
% The problem: a layered column with, per layer, an axis theta, a contrast
% dlam and a scattering ratio r_db. The observables at each depth z and
% synthetic azimuth psi are
%   dP_hh, dP_hv  co- and cross-polarized power anomalies [dB]
%   phi           HH-VV coherence phase [rad]
%   Cmag          |C|, used to size phi's error (lesson 2)
% The forward model is ptt.fujitaModel (Fujita et al. 2006): products of
% 2x2 rotation and transmission matrices down the column and back.

clear; close all; rng(7);
tutorial_setup();

%% 1. A synthetic column (the setup of test_fabric_gls)
fc  = 750e6;
z   = (10:10:1400)';                 % 140 depths
az  = deg2rad(0:10:170);             % 18 azimuths (named az: psi() is MATLAB's trigamma, used below)
nL  = 8;
top = (0:nL-1)' * (1400/nL);         % layer tops [m]
theta_true = deg2rad(35) * ones(nL,1);
dlam_true  = 0.03 + 0.04*(top/1400);
rdb_true   = zeros(nL,1);

lay = struct('top_m', num2cell(top), 'dlam', num2cell(dlam_true), ...
  'theta', num2cell(theta_true), 'r_db', num2cell(rdb_true), ...
  'gx_db', num2cell(zeros(nL,1)));
fm = ptt.fujitaModel(lay, z, az, struct('fc', fc, 'win_m', 30));

figure('Name', 'Lesson 7a', 'Position', [100 100 1200 400]);
flds = {'dP_hh', 'dP_hv', 'phi', 'Cmag'};
ttl  = {'dP_{hh} [dB]', 'dP_{hv} [dB]', '\phi [rad]', '|C|'};
for q = 1:4
  subplot(1,4,q);
  imagesc(rad2deg(az), z, fm.(flds{q})); colorbar;
  if q <= 2, clim([-15 5]); end
  xlabel('\psi [deg]'); if q == 1, ylabel('depth [m]'); end
  title(ttl{q});
end
% The patterns repeat with depth as the phase accumulates; their position
% in psi is set by the fabric axis (lesson 6).

%% 2. Noise and the data covariance (H_pack) - lessons 2, 3
% The standard deviation of a dB power from N looks is known exactly:
%     s_db = (10/log(10)) * sqrt(psi(1, nlook));      psi(1,x) = trigamma
% The phase sigma is lesson 2's CRB.
nlook = 20;
s_db = (10/log(10)) * sqrt(psi(1, nlook))         % ~1 dB for 20 looks
c = min(max(fm.Cmag, 1e-3), 0.999);
sd_phi = sqrt((1 - c.^2) ./ (2*nlook*c.^2));

obs = struct('psi', az, 'Cmag', fm.Cmag, ...
  'dP_hh', fm.dP_hh + s_db*randn(size(fm.dP_hh)), ...
  'dP_hv', fm.dP_hv + s_db*randn(size(fm.dP_hv)), ...
  'phi',   fm.phi   + sd_phi.*randn(size(fm.phi)));

% H_pack stacks every observable into one column d with one sigma each;
% C_d is diagonal (Cdi = 1./sd.^2). Correlation between samples is handled
% by inflation (lesson 3, option c) through n_looks, n_indep_z, n_indep_psi.
numel(z) * numel(az) * 3                         % length of d

%% 3. Unknowns and the prior - lesson 5
% m stacks three profiles on nL layers:
%     ip = struct('th', 1:nL, 'dl', nL+(1:nL), 'rd', 2*nL+(1:nL));
% 24 unknowns here. The prior is lesson 5's exponential covariance per
% profile, combined with blkdiag. Defaults, [mean sigma corr_len_m]:
%     theta_prior [0 pi/2 1e4]      dlam_prior [0.05 0.08 300]      rdb_prior [0 4 300]
z_mid = top + 1400/nL/2;
gp = @(zc, s, L) s^2 * exp(-abs(zc - zc') / L);
Cm = blkdiag(gp(z_mid, pi/2, 1e4), gp(z_mid, 0.08, 300), gp(z_mid, 4, 300));

figure('Name', 'Lesson 7b', 'Position', [100 100 460 400]);
imagesc(Cm ./ sqrt(diag(Cm) * diag(Cm)')); axis image; colorbar;
title('Prior correlation matrix');
xlabel('parameter'); ylabel('parameter');
% The theta block is ~1 throughout: a 10 km correlation length holds one
% axis for the column, as a prior only.

%% 4. The residual (H_resid) - lesson 6
%     r = d - g;
%     r(is_ph) = angle(exp(1i * r(is_ph)));   % wrap the phase entries
% Power anomalies have -Inf nulls on the fabric axes, so data and model
% are floored at clip_db (-30 dB) before subtracting.

%% 5. The loss - lesson 5
%     loss = sum(Cdi .* r.^2) + (m-m_p)' * Cmi * (m-m_p)

%% 6. Multi-start - lesson 6
% Score a 6 x 6 grid of constant-theta, constant-dlam starts without
% iterating, then solve from the best n_start.

%% 7. Gauss-Newton with a prior (H_solve) - lessons 5, 6
%     G  = H_jac(fwd, m, ...)                      % finite-difference Jacobian
%     A  = G'*(Cdi .* G) + Cmi;
%     b  = G'*(Cdi .* r) - Cmi*(m - m_p);
%     dm = A \ b;
% then halve the step until the loss drops. H_jac steps each parameter by
% 1e-3 of its prior sigma, not of its value: r_db starts at 0.

%% 8. Run it
opts = struct('fc', fc, 'n_layer', nL, 'n_looks', nlook, 'max_iter', 30, 'n_start', 3);
tic; out = ptt.fabricGLS(obs, z, opts); t_run = toc;
fprintf('fabricGLS took %.1f s, converged = %d after %d steps\n', t_run, out.converged, out.n_iter);

% The posterior (end of fabricGLS) - lessons 1, 4, 5:
%     C_M  = H_inv(A)
%     Rres = C_M * (G' * (Cdi .* G))      resolution matrix
%     dof  = numel(d_obs) - sum(res)
%     sig  = sig * sqrt(max(chi2_dof, 1)) widen, never shrink
fprintf('chi2/dof = %.2f\n', out.chi2_dof);
fprintf('resolution, theta layers: %s\n', sprintf('%.2f ', out.resolution(out.ip.th)));
fprintf('resolution, dlam layers:  %s\n', sprintf('%.2f ', out.resolution(out.ip.dl)));

% theta is an axis (mod 180): fold onto the truth's branch before comparing.
fold = @(a, ref) ref + angle(exp(2i*(a - ref)))/2;
th_z = fold(out.theta_z, deg2rad(35));
dl_true_z = interp1(top, dlam_true, z, 'previous', 'extrap');

figure('Name', 'Lesson 7c', 'Position', [100 100 1000 420]);
subplot(1,3,1);
band(rad2deg(th_z), rad2deg(out.sigma_theta_z), z); hold on;
plot(rad2deg(th_z), z, 'LineWidth', 1.5); xline(35, 'k', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse'); xlabel('\theta [deg]'); ylabel('depth [m]');
title('\theta(z) \pm2\sigma'); grid on; xlim([25 45]);
subplot(1,3,2);
band(out.dlam_z, out.sigma_dlam_z, z); hold on;
plot(out.dlam_z, z, 'LineWidth', 1.5); plot(dl_true_z, z, 'k', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse'); xlabel('\Delta\lambda'); title('\Delta\lambda(z) \pm2\sigma'); grid on;
subplot(1,3,3);
barh(out.z_mid, [out.resolution(out.ip.th), out.resolution(out.ip.dl), out.resolution(out.ip.rd)]);
set(gca, 'YDir', 'reverse'); xlabel('diag(R)'); xlim([0 1.05]);
legend('\theta', '\Delta\lambda', 'r_{dB}', 'Location', 'southoutside', 'Orientation', 'horizontal');
title('Resolution per layer'); grid on;
% dlam_z is interpolated linearly between layer centres, hence the line
% through a stepped truth. The +-2 sigma bands are thinner than the lines:
% 7560 data at 20 looks constrain a model-generated column tightly.

%% 9. Error-bar check (test_fabric_gls verdict 3)
% 16 noise realizations with the repository test's seeds; compare the
% scatter at mid-column with the reported sigma. About 5 s each.
R = 16;
th_hat = nan(R,1); dl_hat = nan(R,1); s_th = nan(R,1); s_dl = nan(R,1); chi2 = nan(R,1);
mid = round(numel(z)/2);
fprintf('\n  k   dlam       sigma      theta [deg]  sigma    chi2/dof\n');
for k = 1:R
  rng(100 + k);
  ok = struct('psi', az, 'Cmag', fm.Cmag, ...
    'dP_hh', fm.dP_hh + s_db*randn(size(fm.dP_hh)), ...
    'dP_hv', fm.dP_hv + s_db*randn(size(fm.dP_hv)), ...
    'phi',   fm.phi   + sd_phi.*randn(size(fm.phi)));
  ok = ptt.fabricGLS(ok, z, opts);
  th_hat(k) = fold(ok.theta_z(mid), deg2rad(35));
  dl_hat(k) = ok.dlam_z(mid);
  s_th(k) = ok.sigma_theta_z(mid);
  s_dl(k) = ok.sigma_dlam_z(mid);
  chi2(k) = ok.chi2_dof;
  fprintf(' %2d   %.6f   %.6f   %8.4f     %.4f   %6.2f\n', k, dl_hat(k), s_dl(k), ...
    rad2deg(th_hat(k)), rad2deg(s_th(k)), chi2(k));
end

% A chi2/dof far above 1 marks a fit in a false minimum (lesson 6); its
% sigma has already been widened by sqrt(chi2) (lesson 4).
good = chi2 < 2;
ratio = @(sc, est) median(sc) / std(est);
fprintf('\n                      all %2d fits      the %2d with chi2/dof < 2\n', R, nnz(good));
fprintf('theta  claimed/measured   %.2f               %.2f\n', ...
  ratio(s_th, th_hat), ratio(s_th(good), th_hat(good)));
fprintf('dlam   claimed/measured   %.2f               %.2f\n', ...
  ratio(s_dl, dl_hat), ratio(s_dl(good), dl_hat(good)));

% Two effects: (a) on converged fits the theta sigma is about 2x too small,
% the linearized posterior of a 24-parameter nonlinear problem; (b) one
% false minimum inflates the scatter of everything. On the good fits the
% dlam sigma is conservative. The repository test reports the all-fits
% column and documents both sigmas as lower bounds. On real data, the
% per-fit chi2/dof is what exposes case (b). For an interval, resample
% (ptt.quadpolJackknife).

%% TAKEAWAYS: where each lesson appears in fabricGLS
%   lesson 1  normal equations, C_M, Monte Carlo       A \ b, C_M, verdict 3
%   lesson 2  1/variance weights, CRB phase sigma      H_pack, Cdi
%   lesson 3  correlated samples, variance inflation   n_looks, n_indep_*
%   lesson 4  chi2/dof, widen only                     chi2_dof, sqrt(max(chi2,1))
%   lesson 5  prior, posterior, resolution, abstain    Cm, Cmi, resolution
%   lesson 6  Jacobian, Gauss-Newton, multi-start,     H_jac, H_solve, start grid,
%             wrapped residuals, linearized sigma      H_resid, "lower bounds"
%
%% TRY THIS
% 1. Rotate the axis with depth (verdict 4):
%      theta_true = deg2rad(35) + deg2rad(50)*(top/max(top));
%    Rerun with the default theta prior and with theta_prior = [0 pi/2 300].
%    Does the long correlation length stop the data following the rotation?
% 2. Drop the phase: opts.use = {'dP_hh', 'dP_hv'}. Which resolutions fall?
% 3. Set opts.n_looks = 80 on data made with 20. Predict chi2/dof first.

%% ------------------------------------------------------------------------
function band(x, s, z)
% Shaded +-2 sigma band around a depth profile.
fill([x - 2*s; flipud(x + 2*s)], [z; flipud(z)], [0.8 0.87 0.95], 'EdgeColor', 'none');
end
