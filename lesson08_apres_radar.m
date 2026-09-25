%% Lesson 8 - The ApRES measurement
%
% ApRES (Autonomous phase-sensitive Radio-Echo Sounder) is a small
% ground-based ice-penetrating radar (Brennan et al. 2014). This lesson
% covers how it forms a range profile, how the GHOST 2023/24 polarimetric
% sites were measured, and how these data differ from the Ridge A survey
% that fabric_anisotropy was built for.

clear; close all;
tutorial_setup();
root = fullfile(fileparts(which('tutorial_setup')), 'data', 'GHOST24_Polarimetric_pRES_OZ');
site = fullfile(root, 'PpRES_20240551_003');

%% 1. FMCW ranging
% ApRES transmits a chirp: a tone whose frequency rises linearly from
% 200 to 400 MHz over one second. An echo from range R returns after
% tau = 2R/ci, so while it arrives the transmitter has moved on by
% (dF/dt) * tau in frequency. Mixing echo with transmitter gives a beat
% at that frequency difference, proportional to range. The recorded
% signal is the sum of one beat tone per reflector, and an FFT separates
% them by frequency, i.e. by range.
d = dir(fullfile(site, '*_HH_*.dat'));
B = apres.loadBurst(fullfile(site, d(1).name));
fprintf('%s: %d chirps of %.1f s, %.0f-%.0f MHz\n', d(1).name, size(B.vif,1), B.T, B.f0/1e6, B.f1/1e6);

ci = B.ci;                               % wave speed in ice [m/s]
R_ex = 1000;                             % example reflector [m]
tau = 2*R_ex/ci;
f_beat = B.B/B.T * tau;
fprintf('reflector at %d m: delay %.1f us, beat %.0f Hz\n', R_ex, tau*1e6, f_beat);

figure('Name', 'Lesson 8a', 'Position', [100 100 1100 330]);
subplot(1,3,1);
tt = linspace(0, 60e-6, 200);
plot(tt*1e6, B.B*tt/B.T/1e3, 'LineWidth', 1.5); hold on;
plot(tt*1e6, B.B*(tt - tau)/B.T/1e3, '--', 'LineWidth', 1.5);
xlim([0 60]); ylim([-3 12]);
xlabel('time after chirp start [\mus]'); ylabel('frequency above 200 MHz [kHz]'); grid on;
legend('transmitted', sprintf('echo from %d m', R_ex), 'Location', 'northwest');
title('Chirp and delayed echo, first 60 \mus');
% The vertical gap between the lines is the beat frequency, constant for
% the rest of the chirp.

subplot(1,3,2);
t = (0:B.samples_per_chirp-1) / B.fs;
plot(t, B.vif(1, :)); xlabel('time in the chirp [s]'); ylabel('voltage [V]');
title('Recorded beat signal, one chirp'); grid on; xlim([0 1]);

%% 2. From chirps to a range profile
% 100 chirps are averaged (noise falls as 1/sqrt(100)), windowed, and
% Fourier transformed (apres.rangeProcess). The phase of each range bin is
% referenced so that it measures the reflector's position to a fraction
% of a wavelength: the "phase-sensitive" part of the name, used for
% time-lapse strain and melt. Here we use the complex profile itself.
P = apres.rangeProcess(B);
subplot(1,3,3);
plot(20*log10(abs(P.spec)), P.range); set(gca, 'YDir', 'reverse');
xlabel('power [dB]'); ylabel('range in ice [m]'); grid on;
title(sprintf('Range profile, %d chirps averaged', P.n_chirps));
% Range bins are 0.21 m apart. The bright return near 2050 m is the bed.

%% 3. A polarimetric site: four shots
% ApRES has one transmit and one receive antenna, each a dipole that
% responds to one polarization. For a polarimetric measurement the
% operator records four shots, turning the antennas between them
% (GHOST24_PpRES_Log_OZ.pdf):
%
%     shot   Tx antenna   Rx antenna
%     HH     H            H
%     HV     H            V
%     VH     V            H
%     VV     V            V
%
% The Tx-Rx line ran true south to north with the antennas ~4 m apart, and
% the H antennas lay along that line, so H points north. Each shot is a
% separate acquisition a few minutes apart. Three consequences:
%   - the antennas move by centimetres between shots, so each shot has its
%     own phase offset (lesson 9 removes it);
%   - the four range axes differ by a fraction of a bin (lesson 9 aligns them);
%   - the only averaging available is over depth: one site has no
%     along-track dimension.
Q = apres.loadQuadpolSite(site, struct('coregister', false));
figure('Name', 'Lesson 8b', 'Position', [100 100 900 400]);
k = ones(47, 1) / 47;                                 % ~10 m running mean
pols = {'hh', 'hv', 'vh', 'vv'};
for p = 1:4
  plot(10*log10(conv(abs(Q.(pols{p})).^2, k, 'same')), Q.z, 'LineWidth', 1); hold on;
end
set(gca, 'YDir', 'reverse'); ylim([0 2200]); grid on;
xlabel('power [dB], 10 m average'); ylabel('depth [m]');
legend(upper(pols), 'Location', 'northwest'); title('Power vs depth, four shots at site 003');
% The cross-polarized shots (HV, VH) are weaker than the co-polarized
% ones, and nearly identical to each other, as reciprocity requires.

%% 4. How this differs from the Ridge A survey
% fabric_anisotropy was written for the CReSIS ground-based accumulation
% radar towed along survey lines at Ridge A. The estimators are the same;
% the data are not.
%
%                      ApRES (GHOST sites)        Ridge A (CReSIS accum radar)
%   centre frequency   300 MHz                    750 MHz
%   phase per unit     0.12 rad/m                 0.30 rad/m
%     dlam
%   acquisition        stationary; four shots     moving; four channels recorded
%                      minutes apart              together, every trace
%   what varies        depth only                 depth and position along track
%   averaging (looks)  depth window only          along-track blocks (125 traces)
%                                                 plus depth
%   channel alignment  measured per site from     coregistered by the OPR toolbox
%                      the profiles               across the image
%   instrument         antenna phase offsets      a cross-pol pedestal in the
%     artefact         between shots              antenna frame
%   azimuth reference  antenna heading from the   vehicle heading, changes along
%                      field log (H = north)      the line
%   product            one column per site        theta and dlam along the line,
%                                                 mapped
%   inversion used     the windowed estimator     ptt.fabricGLS / quadpolFabricLS
%                      (project)
%
% Two of these matter most for the estimator. At 300 MHz the phase turns
% 2.5x more slowly per metre, so a window sees less of the pattern and
% weak fabric (dlam < ~0.07) cannot be resolved in 60 m. And with no
% along-track looks, the number of independent samples per window is
% small (about 10 per 10 m), which is why the error bars must be checked
% by resampling rather than trusted from the formula.

%% TRY THIS
% 1. The bed at site 003 is near 2050 m. What beat frequency is that?
%    (Section 1's formula.) What sets the maximum range ApRES can record?
% 2. Compare the HV and VH profiles with conv(hv.*conj(vh)). How similar
%    are they, and why should they be?
