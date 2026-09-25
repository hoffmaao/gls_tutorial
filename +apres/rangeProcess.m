function P = rangeProcess(B, opts)
%RANGEPROCESS Turn a burst's deramped chirps into a complex range profile.
%
% P = apres.rangeProcess(B, opts)
%
% Port of fmcw_range_vel (Brennan et al. 2013, via a Python translation
% it agrees with to 1e-7), on the chirp-averaged burst.
%
% HOW FMCW RANGING WORKS, in one paragraph. The radar sweeps its frequency
% linearly from f0 to f1 over T seconds and mixes the echo with what it is
% transmitting at that moment. An echo from range R arrives tau = 2R/ci
% late, so the two differ by a CONSTANT frequency K*tau/(2 pi): the "beat"
% or "deramped" signal vif is a sum of sine waves, one per reflector, with
% frequency proportional to range. An FFT separates sine waves by
% frequency, so the FFT of vif IS the range profile. Its phase is then
% referenced (phiref, Brennan eq. 17) so it measures the fine position of
% each reflector - the "phase-sensitive" in pRES/ApRES.
%
% The FFT is also just a matrix times a vector (lesson 0): every range bin
% is one row of the DFT matrix dotted with the chirp. Nothing in this
% function is nonlinear except the final phase-to-range conversion.
%
% opts
%   .pad        (2)     zero-padding factor
%   .max_range  (3000)  metres of profile to keep
%   .att_index  (1)     which attenuator setting to use
%
% P fields
%   range   [Nr x 1] coarse range bin centres (m)
%   spec    [Nr x 1] phase-referenced complex profile (the one to use)
%   spec_raw[Nr x 1] before the phase reference
%   n_chirps number of chirps averaged

if nargin < 2, opts = struct(); end
pad  = opt(opts, 'pad', 2);
rmax = opt(opts, 'max_range', 3000);
ai   = opt(opts, 'att_index', 1);

% --- one attenuator setting (fmcw_burst_split_by_att), then the mean chirp
atts = unique(B.chirp_att, 'stable');
ai = min(max(ai, 1), numel(atts));
keep = B.chirp_att == atts(ai);
vif = mean(B.vif(keep, :), 1);            % averaging chirps: noise down by sqrt(n)

% --- coarse range axis (func_depth_conversion), cropped at rmax
N  = B.samples_per_chirp;
nf = round(pad*N/2 - 0.5);
Rc = (0:nf-1)' * B.ci / (2 * B.B * pad);
Rc = Rc(Rc <= rmax);
nr = numel(Rc);
n  = (0:nr-1)';

% --- window, zero-pad, rotate so the phase is measured at t = T/2, FFT
k   = (0:N-1)';
win = 0.42 - 0.5*cos(2*pi*k/(N-1)) + 0.08*cos(4*pi*k/(N-1));   % numpy blackman
x   = vif(:) - mean(vif);
xw  = win .* x;
nfft = pad * N;
vp  = zeros(nfft, 1);
vp(1:N) = xw;
if mod(N, 2) == 1, xn = (N-1)/2; else, xn = N/2; end
vp  = circshift(vp, -xn);
F   = sqrt(2*pad) / nfft * fft(vp) / sqrt(mean(win.^2));
s   = F(1:nr);

% --- phase reference (Brennan eq. 17)
phiref = 2*pi*B.fc*n/(B.B*pad) - B.K*n.^2/(2*B.B^2*pad^2);
P = struct('range', Rc, 'spec', exp(-1i*phiref) .* s, 'spec_raw', s, ...
  'n_chirps', nnz(keep), 'att', atts(ai), 'timestamp', B.timestamp, ...
  'fc', B.fc, 'B', B.B, 'ci', B.ci, 'lambdac', B.lambdac, 'file', B.file);
end

function v = opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
