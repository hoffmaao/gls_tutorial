%% Lesson 9 - From ApRES files to the coherence field
%
% The project's estimator fits the HH-VV coherence as a function of depth
% and antenna azimuth. This lesson builds that field for one GHOST site,
% one function at a time. Each function in +apres/ explains its own step
% in its header.

clear; close all;
tutorial_setup();
root = fullfile(fileparts(which('tutorial_setup')), 'data', 'GHOST24_Polarimetric_pRES_OZ');
site = fullfile(root, 'PpRES_20240551_003');

%% 1. The four shots, aligned (apres.loadQuadpolSite)
% Each shot is loaded and range processed (lesson 8). The four range
% axes differ by a fraction of a bin because the antennas were moved; each
% channel's offset is measured against HH from the power envelope and
% removed by a sub-bin shift. Misaligned channels decorrelate.
Q = apres.loadQuadpolSite(site);
fprintf('range offsets vs HH [bins]: HV %.2f, VH %.2f, VV %.2f\n', ...
  Q.offset_bins.hv, Q.offset_bins.vh, Q.offset_bins.vv);

%% 2. Antenna phase offsets (apres.calibratePhase)
% Moving an antenna by centimetres changes its phase by tens of degrees.
% In the top 40 m there is no fabric, so HH and VV should be in phase. Two
% measurements fix the transmit and receive antenna phases: reciprocity
% (HV must equal VH) and the near-surface HH-VV phase. One 180-degree
% ambiguity remains; see the function header.
k = ones(47, 1) / 47;                                 % ~10 m running mean
hhvv = @(Q) conv(Q.hh .* conj(Q.vv), k, 'same');
before = angle(hhvv(Q));
Q = apres.calibratePhase(Q);
after = angle(hhvv(Q));
near = Q.z > 10 & Q.z < 40;
fprintf('HH-VV phase at 10-40 m: %.0f deg before, %.0f deg after\n', ...
  rad2deg(angle(mean(exp(1i*before(near))))), rad2deg(angle(mean(exp(1i*after(near))))));
fprintf('calibration: a = %.0f deg, b = %.0f deg (reciprocity coherence %.2f)\n', ...
  rad2deg(Q.phase_cal.a), rad2deg(Q.phase_cal.b), Q.phase_cal.coh_reciprocity);

%% 3. Synthesized azimuths and the coherence field (apres.coherenceField)
% Rotating the antennas by psi is a linear combination of the four
% channels (ptt.rotatePolarization), so the response at any azimuth
% follows from one measurement. Coherence is averaged over 10 m windows;
% the number of independent samples per window is measured from the
% speckle autocorrelation (lesson 3). psi is clockwise from the H antenna,
% which points north at these sites.
F = apres.coherenceField(Q, struct('z_range', [20 1200]));
fprintf('%d depths x %d azimuths, %.1f independent looks per 10 m window\n', ...
  numel(F.z), numel(F.psi), F.n_looks);

figure('Name', 'Lesson 9', 'Position', [100 100 1000 420]);
subplot(1,3,1);
imagesc(rad2deg(F.psi), F.z, abs(F.C)); colorbar; clim([0 1]);
xlabel('antenna azimuth \psi [deg from north]'); ylabel('depth [m]'); title('|C|');
subplot(1,3,2);
imagesc(rad2deg(F.psi), F.z, angle(F.C)); colorbar;
xlabel('\psi [deg from north]'); title('arg C [rad]');
subplot(1,3,3);
imagesc(rad2deg(F.psi), F.z, F.sigma); colorbar;
xlabel('\psi [deg from north]'); title('\sigma per real/imag part');
% The estimator fits the phase panel: the position of the pattern in psi
% gives theta0, its change with depth gives ddelta and so dlam (lesson 6).

%% 4. How precise is each row? The data covariance
% All 18 azimuths in a row are made from the same four channels, so their
% errors are shared. apres.coherenceCov computes the covariance of
% [Re C; Im C] for one row from its 4x4 channel moments M and the number
% of looks N. Check it the lesson-1 way: simulate many N-look
% measurements of one row and compare.
i = find(F.z > 600, 1);                              % one row near 600 m
S = struct('hh', Q.hh, 'vv', Q.vv, 'hv', Q.hv, 'vh', Q.vh);
Mall = ptt.quadpolMoments(S, [47 1]);                % 10 m of 0.21 m bins
[~, iz] = min(abs(Q.z - F.z(i)));
M = reshape(Mall(iz, :, :), 4, 4);
psi_s = (0:10:170) * pi/180;                         % synthesis-sense azimuths
N = round(F.n_looks);
Cd = apres.coherenceCov(M, N, psi_s);

R = 5000; X = zeros(numel(Cd(:,1)), R);
L = chol(M + 1e-12*trace(M)*eye(4), 'lower');
for r = 1:R
  x = L * (randn(4, N) + 1i*randn(4, N)) / sqrt(2);  % N looks of the 4 channels
  [~, Cr] = apres.coherenceCov((x*x') / N, N, psi_s);
  X(:, r) = [real(Cr); imag(Cr)];
end
ratio = sqrt(diag(Cd)) ./ std(X, 0, 2);
fprintf('row at %.0f m, %d looks: predicted/simulated sigma, median %.2f (range %.2f-%.2f)\n', ...
  F.z(i), N, median(ratio), min(ratio), max(ratio));
ev = sort(eig(Cd), 'descend');
fprintf('independent numbers in the row: %d of %d\n', nnz(ev > 1e-9*ev(1)), numel(ev));

figure('Name', 'Lesson 9b', 'Position', [100 100 900 380]);
subplot(1,2,1);
imagesc(Cd ./ sqrt(diag(Cd) * diag(Cd)')); axis image; colorbar; clim([-1 1]);
xlabel('entry of [Re C; Im C]'); ylabel('entry of [Re C; Im C]');
title(sprintf('Correlation of one row''s errors, %.0f m', F.z(i)));
subplot(1,2,2);
semilogy(ev / ev(1), 'o-', 'LineWidth', 1.2); grid on;
xlabel('eigenvalue index'); ylabel('eigenvalue / largest');
title('Eigenvalues of the row covariance');
% Most eigenvalues are zero to machine precision: 36 numbers per row, about
% 7 of them independent. A fit that treats all 36 as independent counts the
% same information several times over and reports error bars that are too
% small (lesson 3). apres.coherenceField keeps the independent directions
% and hands the estimator a whitening matrix for each row (F.Wh).
%
% Two more fields matter for the fit. F.z_eff: each row averages 10 m of
% ice weighted by random speckle brightness, so its coherence belongs to
% the brightness-weighted depth, not the row centre. And the smallest
% eigenvalues kept are limited to 5% of the largest (cond_floor), because
% directions the covariance calls near-perfect amplify any model error.

%% 5. A synthetic site (apres.syntheticSite)
% The same chain runs on a site with a chosen fabric, the known answer
% the estimator is tested against first. syntheticSite shifts VV by 0.6
% bins and co-registers it like a real site.
[Qs, truth] = apres.syntheticSite(struct('theta', deg2rad(30) * ones(5, 1)));
fprintf('synthetic VV range offset: %.2f bins injected, %.2f measured\n', ...
  Qs.true_offset_bins.vv, Qs.offset_bins.vv);
Qs = apres.calibratePhase(Qs);
Fs = apres.coherenceField(Qs, struct('z_range', [20 1200]));
fprintf('synthetic: axis 30 deg, dlam %s by layer\n', mat2str(truth.dlam', 3));

%% TRY THIS
% 1. Repeat for sites 001, 002 and 004. Which has the lowest reciprocity
%    coherence, and what does that mean for its calibration?
% 2. Skip calibratePhase. What happens to the phase panel near the surface?
