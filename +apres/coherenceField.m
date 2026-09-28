function F = coherenceField(Q, opts)
%COHERENCEFIELD The HH-VV coherence at every depth and antenna azimuth.
%
% F = apres.coherenceField(Q, opts)
%
% The data the project's estimator fits. Q is a calibrated site
% (apres.loadQuadpolSite, then apres.calibratePhase) or a synthetic one
% (apres.syntheticSite, then apres.calibratePhase).
%
% This is apres.observables with one change of convention. observables
% synthesizes in the sense ptt.ershadiFabric and ptt.fabricGLS use, where
% a fabric axis at azimuth theta shows its pattern at psi = -theta. Here
% psi is flipped so the pattern sits at psi = +theta: fit
% mu = cos 2(psi - theta0) and theta0 IS the axis azimuth, clockwise from
% the H antenna. On the GHOST sites H points to true north (field log),
% so theta0 is a compass bearing.
%
% opts: as apres.observables (defaults win_m 10, psi_step_deg 10), plus
%       .cond_floor (0.05) smallest variance kept, relative to the largest
%
% F fields: z [Nz x 1] row centres, z_eff [Nz x 1] power-weighted depth of
%   each row (the depth its coherence belongs to), psi [1 x Np] (rad,
%   0..pi), C [Nz x Np] complex coherence, Cd [2Np x 2Np x Nz] covariance of [Re C; Im C] per row,
%   Wh [2Np x 2Np x Nz] whitening matrix per row (zero rows beyond rank),
%   rank [Nz x 1] independent numbers per row, sigma [Nz x Np] typical
%   error per real/imaginary part, n_looks, snr_db [Nz x 1], fc, win_m,
%   dz (range-bin spacing, m)

if nargin < 2, opts = struct(); end
if ~isfield(opts, 'win_m'), opts.win_m = 10; end
O = apres.observables(Q, opts);
[psi, order] = sort(mod(-O.psi, pi));
Np = numel(psi);
idx = [order, Np + order];
Cd = O.Cd(idx, idx, :);

% WHITENING (lesson 3). Each row's covariance is singular: 18 azimuths made
% from four channels carry only about 7 independent numbers. Keep the
% directions the data measure and scale each to unit variance, so that
% r = Wh * [Re e; Im e] has independent, unit-variance entries.
% CAP. The exact covariance claims some directions are known almost
% perfectly, where any small model error then dominates the fit (tested
% near phase nodes on synthetic sites). No direction's variance is allowed
% below cond_floor (default 0.05) times the row's largest; this errs on
% the cautious side.
Nz = numel(O.z);
Wh = zeros(2*Np, 2*Np, Nz);
rank = zeros(Nz, 1);
fl = 0.05;
if isfield(opts, 'cond_floor'), fl = opts.cond_floor; end
for i = 1:Nz
  [U, L] = eig(Cd(:, :, i), 'vector');
  k = L > 1e-9 * max(L);
  rank(i) = nnz(k);
  L = max(L, fl * max(L));
  Wh(1:rank(i), :, i) = diag(1 ./ sqrt(L(k))) * U(:, k).';
end

F = struct('z', O.z, 'z_eff', O.z_eff, 'psi', psi, 'C', O.C(:, order), 'Cd', Cd, 'Wh', Wh, ...
  'rank', rank, 'sigma', O.sigma(:, order), 'n_looks', O.n_looks, ...
  'snr_db', O.snr_db, 'fc', O.fc, 'win_m', O.win_m, 'dz', O.dz);
end
