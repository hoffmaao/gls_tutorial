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
% opts: as apres.observables; defaults win_m 10, psi_step_deg 10.
%
% F fields: z [Nz x 1], psi [1 x Np] (rad, 0..pi), C [Nz x Np] complex
%   coherence, sigma [Nz x Np] noise per real/imaginary part, n_looks,
%   snr_db [Nz x 1], fc, win_m

if nargin < 2, opts = struct(); end
if ~isfield(opts, 'win_m'), opts.win_m = 10; end
O = apres.observables(Q, opts);
[psi, order] = sort(mod(-O.psi, pi));
F = struct('z', O.z, 'psi', psi, 'C', O.C(:, order), ...
  'sigma', O.sigma(:, order), 'n_looks', O.n_looks, 'snr_db', O.snr_db, ...
  'fc', O.fc, 'win_m', O.win_m);
end
