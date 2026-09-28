function O = observables(Q, opts)
%OBSERVABLES Quad-pol ApRES profiles -> the observables ptt.fabricGLS fits.
%
% O = apres.observables(Q, opts)
%
% Q is what apres.loadQuadpolSite returns (hh, hv, vh, vv, z, fc). The
% output O can be handed straight to ptt.fabricGLS(O, O.z, ...).
%
% THE STEPS (lesson 9 walks through each one):
%   1. MOMENTS. Every observable is quadratic in the scattering matrix
%      (a power or a correlation), so first form the 16 averages
%      <S_k conj(S_l)> over a depth window of win_m (ptt.quadpolMoments
%      with range looks). Averaging over depth is what turns one noisy
%      speckle sample into a measurement with N looks - an ApRES site has
%      no along-track traces to average over, so depth is all we have.
%   2. SYNTHESIS. Rotating the antennas by psi is linear in S, so each
%      synthetic power and correlation is a fixed trig combination of the
%      moments. Same weights and same rotation sense as ptt.ershadiFabric
%      (its E1), which is the sense ptt.fujitaModel and ptt.fabricGLS use.
%   3. OBSERVABLES. dP_hh, dP_hv: power anomalies in dB against the
%      azimuthal mean amplitude (Ershadi eq. 12). C: HH-VV coherence from
%      the same moments, conjugated for deramped FMCW data (their E2).
%   4. LOOKS. The number of INDEPENDENT samples in the window is not
%      win_m / dz: range bins are oversampled (zero padding) and smeared
%      by the window function, so neighbours are correlated (lesson 3).
%      It is measured from the data's own speckle autocorrelation:
%          n_looks = nr / sum_k |rho_k|^2
%      where rho_k is the correlation of the complex profile with itself
%      shifted by k bins. This is the n_looks fabricGLS needs.
%   5. DECIMATION. Rows closer than the window share samples, so their
%      errors are correlated. Keeping one row per win_m makes the rows
%      independent, so C_d is block diagonal: one full block per row
%      (lesson 3).
%
% opts
%   .win_m        (20)   depth window for the moments [m]
%   .psi_step_deg (10)   synthetic azimuth step [deg]
%   .z_range      ([])   [zmin zmax] depths to keep; default all
%   .deramped     (true) FMCW data are deramped: conjugate C
%
% O fields
%   z [Nz x 1], psi [1 x Np], dP_hh, dP_hv, phi, Cmag [Nz x Np],
%   C [Nz x Np] complex HH-VV coherence, Cd [2Np x 2Np x Nz] covariance of
%   [Re C; Im C] for each row, sigma [Nz x Np] typical error per real or
%   imaginary part (for plots),
%   n_looks, win_m, fc, snr_db [Nz x 1] (co-pol power over the noise floor),
%   rho (the speckle autocorrelation used for n_looks)

if nargin < 2, opts = struct(); end
win_m = opt(opts, 'win_m', 20);
step  = opt(opts, 'psi_step_deg', 10);
zr    = opt(opts, 'z_range', []);
deramped = opt(opts, 'deramped', true);

z = Q.z(:);
dz = z(2) - z(1);
nr = 2*floor(win_m / dz / 2) + 1;               % odd number of bins in the window

% --- 1. moments over the depth window
S = struct('hh', Q.hh(:), 'vv', Q.vv(:), 'hv', Q.hv(:), 'vh', Q.vh(:));
M = ptt.quadpolMoments(S, [nr 1]);

% --- 2. synthesis, ptt.ershadiFabric's R S R' sense
psi = (0:step:180-step) * pi/180;
Np = numel(psi); Nt = numel(z);
Phh = zeros(Nt, Np); Pvv = zeros(Nt, Np); Phv = zeros(Nt, Np);
Chv = complex(zeros(Nt, Np));
for j = 1:Np
  c = cos(-psi(j)); s = sin(-psi(j)); cs = c*s;
  w_hh = [c*c, s*s,  cs,  cs];                  % weights on (hh, vv, hv, vh)
  w_vv = [s*s, c*c, -cs, -cs];
  w_hv = [-cs,  cs, c*c, -s*s];
  Phh(:, j) = real(quad(M, w_hh, w_hh));
  Pvv(:, j) = real(quad(M, w_vv, w_vv));
  Phv(:, j) = real(quad(M, w_hv, w_hv));
  Chv(:, j) = quad(M, w_hh, w_vv);
end

% --- 3. observables
A_hh = sqrt(max(Phh, 0));
A_hv = sqrt(max(Phv, 0));
dP_hh = 20*log10(max(A_hh, realmin) ./ max(mean(A_hh, 2), realmin));
dP_hv = 20*log10(max(A_hv, realmin) ./ max(mean(A_hv, 2), realmin));
C = Chv ./ max(sqrt(Phh .* Pvv), realmin);
if deramped, C = conj(C); end

% --- 4. independent looks, from the speckle autocorrelation of HH
[n_looks, rho] = effectiveLooks(S.hh, nr);

% --- signal-to-noise: the deepest 5% of the record is taken as noise
Pco = mean(Phh + Pvv, 2) / 2;
noise = median(Pco(round(0.95*Nt):end));
snr_db = 10*log10(Pco / noise);

% --- 5. one row per window, inside z_range
keep = (1:nr:Nt)';
keep = keep(keep > (nr-1)/2 & keep <= Nt - (nr-1)/2);   % whole windows only
if ~isempty(zr), keep = keep(z(keep) >= zr(1) & z(keep) <= zr(2)); end

% --- effective depth of each row. The moments average the window with
% weights set by the speckle brightness, so a row's coherence belongs to
% the power-weighted mean depth of its window, not to its centre.
Pbin = abs(S.hh).^2 + abs(S.vv).^2 + abs(S.hv).^2 + abs(S.vh).^2;
kr = ones(nr, 1);
z_eff = conv(Pbin .* z, kr, 'same') ./ max(conv(Pbin, kr, 'same'), realmin);

% --- 6. the covariance of the coherence field, row by row
% Every synthesized coherence is a function of the same 4x4 moments, so
% the 18 azimuths in a row share their errors. For N looks of complex
% Gaussian channels the moments' covariance is known exactly, and
% linearizing carries it to [Re C; Im C] (apres.coherenceCov). Rows one
% window apart are independent, so C_d is block diagonal, one block per row.
Ck = C(keep, :);
Nk = numel(keep);
Cd = zeros(2*Np, 2*Np, Nk);
for i = 1:Nk
  Cd(:, :, i) = apres.coherenceCov(reshape(M(keep(i), :, :), 4, 4), n_looks, psi);
end
if deramped                                    % C -> conj(C) flips Im
  Cd(1:Np, Np+1:end, :) = -Cd(1:Np, Np+1:end, :);
  Cd(Np+1:end, 1:Np, :) = -Cd(Np+1:end, 1:Np, :);
end
% typical error of one real or imaginary part, for plotting
sigma = zeros(Nk, Np);
for i = 1:Nk
  dg = diag(Cd(:, :, i));
  sigma(i, :) = sqrt((dg(1:Np) + dg(Np+1:end)).' / 2);
end

O = struct('z', z(keep), 'z_eff', z_eff(keep), 'psi', psi, 'C', Ck, 'sigma', sigma, 'Cd', Cd, ...
  'dP_hh', dP_hh(keep, :), 'dP_hv', dP_hv(keep, :), ...
  'phi', angle(C(keep, :)), 'Cmag', abs(C(keep, :)), ...
  'n_looks', n_looks, 'win_m', win_m, 'fc', Q.fc, ...
  'snr_db', snr_db(keep), 'rho', rho, 'nr', nr, 'dz', dz);
end

% -------------------------------------------------------------------------
function [n, rho] = effectiveLooks(s, nr)
% Number of independent samples in an nr-bin window of a speckle profile.
% For a window average of correlated samples, var(mean) = var/n with
%   n = nr / sum_{k=-(nr-1)}^{nr-1} (1 - |k|/nr) |rho_k|^2
% (lesson 3's triangle). rho is estimated from the whole profile after
% dividing out the slow power trend, so the decay of power with depth does
% not masquerade as correlation.
s = s(:);
p = movmean(abs(s).^2, 5*nr);
u = s ./ sqrt(max(p, realmin));
K = nr - 1;
rho = zeros(K+1, 1);
for k = 0:K
  rho(k+1) = abs(sum(u(1:end-k) .* conj(u(1+k:end)))) / sum(abs(u).^2);
end
w = 1 - (0:K)'/nr;
n = nr / (1 + 2*sum(w(2:end) .* rho(2:end).^2));
end

function v = quad(M, w1, w2)
v = complex(zeros(size(M, 1), 1));
for k = 1:4
  for l = 1:4
    v = v + (w1(k) * w2(l)) * M(:, k, l);
  end
end
end

function v = opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
