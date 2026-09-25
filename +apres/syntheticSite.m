function [Q, truth] = syntheticSite(truth, opts)
%SYNTHETICSITE A fake quad-pol ApRES site whose fabric you know.
%
% [Q, truth] = apres.syntheticSite(truth, opts)
%
% Returns Q in the same form apres.loadQuadpolSite returns real data,
% so everything downstream can be tested on a known answer first.
%
% HOW IT IS MADE.
%   1. ptt.fujitaModel gives the four channel amplitudes (hh, vv, hv = vh)
%      of a layered fabric column, with the antennas at azimuth 0.
%   2. SPECKLE. Real ice returns a random sum of many small reflections in
%      every range bin. Every channel sees the same reflectors, so each
%      depth gets ONE random complex amplitude shared by all four
%      channels, smoothed over a few bins like a real range response.
%   3. POWER DECAY with depth (attenuation and spreading).
%   4. RECEIVER NOISE, independent in each channel.
%   5. INSTRUMENT EFFECTS seen on the GHOST 2023/24 sites: antenna phase
%      offsets between shots (the antennas are re-laid; see
%      apres.calibratePhase for the model) and a sub-bin range offset.
%   6. DERAMPING. Real FMCW data carry the conjugate of the model's
%      coherence (Ershadi et al. 2022; ptt.ershadiFabric note E2), so the
%      channels are conjugated to look like real data.
%   7. CO-REGISTRATION (apres.coregister), exactly as for real data:
%      the offset is measured and removed.
%
% truth (all optional)
%   .top_m   layer tops [m]            default (0:250:1000)'
%   .theta   axis azimuth per layer, radians clockwise from the H antenna
%   .dlam    horizontal eigenvalue difference per layer
%   .r_db    anisotropic scattering ratio per layer (dB), default 0
% opts (all optional)
%   .fc (300e6), .dz (0.21), .z_max (1400)
%   .decay_m (250)        power e-folding depth
%   .snr_surface_db (60)  co-pol power over noise near the surface
%   .antenna_phase ([0.7 -0.4])  [a b]: receive and transmit V-minus-H
%                         antenna phases [rad] (apres.calibratePhase)
%   .range_offset (0.6)   injected VV range offset vs HH [bins]
%   .coregister (true)    measure and remove the range offsets
%   .seed (1)

if nargin < 1, truth = struct(); end
if nargin < 2, opts = struct(); end
top   = opt(truth, 'top_m', (0:250:1000)');
nL    = numel(top);
theta = opt(truth, 'theta', deg2rad(35) * ones(nL, 1));
dlam  = opt(truth, 'dlam', linspace(0.03, 0.12, nL)');
r_db  = opt(truth, 'r_db', zeros(nL, 1));
truth = struct('top_m', top(:), 'theta', theta(:), 'dlam', dlam(:), 'r_db', r_db(:));

fc   = opt(opts, 'fc', 300e6);
dz   = opt(opts, 'dz', 0.21);
zmax = opt(opts, 'z_max', 1400);
rng(opt(opts, 'seed', 1));

z = (dz:dz:zmax)';
lay = struct('top_m', num2cell(truth.top_m), 'dlam', num2cell(truth.dlam), ...
  'theta', num2cell(truth.theta), 'r_db', num2cell(truth.r_db), ...
  'gx_db', num2cell(zeros(nL, 1)));
% eps_perp must match the constant fit_column converts phase rate to dlam
% with (ptt.constants eps_bar = 3.171); fujitaModel's own default is the
% paper's 3.15, which would put a 0.3% scale bias into the truth.
fm = ptt.fujitaModel(lay, z, 0, struct('fc', fc, 'eps_perp', 3.171));

% speckle: one complex Gaussian amplitude per bin, smoothed over ~5 bins
k = [0.2 0.6 1 0.6 0.2]';
a = conv(complex(randn(numel(z)+4, 1), randn(numel(z)+4, 1)), k / norm(k), 'valid');
amp = a .* exp(-z / (2 * opt(opts, 'decay_m', 250)));

sig_n = 10^(-opt(opts, 'snr_surface_db', 60) / 20);
noise = @() sig_n * complex(randn(numel(z), 1), randn(numel(z), 1)) / sqrt(2);

Q = struct();
Q.hh = amp .* fm.s_hh(:, 1) + noise();
Q.vv = amp .* fm.s_vv(:, 1) + noise();
Q.hv = amp .* fm.s_hv(:, 1) + noise();
Q.vh = amp .* fm.s_hv(:, 1) + noise();

% instrument effects (apres.calibratePhase explains the a/b model), then
% the deramp conjugate
ab = opt(opts, 'antenna_phase', [0.7 -0.4]);        % [a b], rad
Q.hv = Q.hv * exp(1i * ab(1));
Q.vh = Q.vh * exp(1i * ab(2));
Q.vv = shiftBins(Q.vv, opt(opts, 'range_offset', 0.6)) * exp(1i * sum(ab));
for p = {'hh', 'vv', 'hv', 'vh'}
  Q.(p{1}) = conj(Q.(p{1}));
end

Q.z = z; Q.fc = fc; Q.site = 'synthetic';
Q.true_offset_bins = struct('hv', 0, 'vh', 0, 'vv', opt(opts, 'range_offset', 0.6));
Q.true_antenna_phase = ab;
if opt(opts, 'coregister', true), Q = apres.coregister(Q, opts); end
end

function y = shiftBins(x, sh)
N = numel(x);
kk = [0:ceil(N/2)-1, -floor(N/2):-1]';
y = ifft(fft(x) .* exp(-2i*pi*kk*sh/N));
end

function v = opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
