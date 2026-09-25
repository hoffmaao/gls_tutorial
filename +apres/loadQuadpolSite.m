function Q = loadQuadpolSite(folder, opts)
%LOADQUADPOLSITE Load the four polarizations of one quad-pol ApRES site.
%
% Q = apres.loadQuadpolSite(folder)
% Q = apres.loadQuadpolSite(folder, opts)
%
% A quad-pol ApRES site is four ordinary ApRES shots taken with the
% antennas turned between them: HH, HV, VH, VV. The polarization is only
% recorded in the FILE NAME (the header cannot know how the antennas were
% laid on the snow). The GHOST 2023/24 files (Zeising) name it, followed
% by the antenna-direction code from the field log:
%   Survey_2024-01-02_083812_HH_RL.dat    (HH | RL, HV | RD, VV | UD, VH | UL)
% so the polarization is taken as the last HH/HV/VH/VV token, optionally
% followed by a two-letter direction code, before the extension, and anything that does not match exactly one file per
% polarization is an error rather than a guess.
%
% Each shot is loaded (apres.loadBurst) and range processed
% (apres.rangeProcess) with the same settings, so all four share one
% range axis.
%
% opts: passed to apres.loadBurst and apres.rangeProcess (burst,
%       samples_per_chirp, er_ice, pad, max_range), plus
%       .att_index   (1) which of the attenuator settings SHARED by all
%                    four shots to use
%       .coregister  (true) sub-sample range alignment to HH, measured
%                    between coreg_zmin (40) and coreg_zmax (1000) m
%
% Q fields
%   hh, hv, vh, vv   [Nr x 1] complex, phase-referenced profiles
%   z                [Nr x 1] depth (m) = range in ice
%   fc               centre frequency (Hz)
%   files            struct of the four source files
%   time             struct of the four shot times
%   site             folder name
%   attenuator       the shared attenuator code used (Attenuator1 + i*AFGain)
%   offset_bins      measured range offset of hv, vh, vv vs hh (bins)

if nargin < 2, opts = struct(); end
d = [dir(fullfile(folder, '*.dat')); dir(fullfile(folder, '*.DAT'))];
d = d(~startsWith({d.name}, '._'));
[~, iu] = unique(lower({d.name}));       % case-insensitive filesystems list twice
d = d(iu);

pols = {'hh', 'hv', 'vh', 'vv'};
Q = struct();
B = cell(1, 4);
for k = 1:numel(pols)
  hit = false(numel(d), 1);
  for j = 1:numel(d)
    [~, stem] = fileparts(d(j).name);
    tok = regexp(stem, '(?<=[-_])(HH|HV|VH|VV)(?=([-_][A-Z]{2})?$)', 'match', 'once', 'ignorecase');
    hit(j) = strcmpi(tok, pols{k});
  end
  if nnz(hit) ~= 1
    error('apres:loadQuadpolSite:pol', ...
      'expected exactly one %s file in %s, found %d', upper(pols{k}), folder, nnz(hit));
  end
  Q.files.(pols{k}) = fullfile(folder, d(hit).name);
  B{k} = apres.loadBurst(Q.files.(pols{k}), opts);
end

% ATTENUATOR SETTINGS MUST MATCH. Rotating the antennas in software mixes
% the four channels (lesson 9), so a gain difference between shots becomes
% a fake fabric signal. Field practice sometimes records one shot with
% an extra setting, so use a setting that ALL four shots recorded; if
% there is none, refuse rather than guess a calibration.
common = B{1}.chirp_att;
for k = 2:4, common = intersect(common, B{k}.chirp_att, 'stable'); end
common = unique(common, 'stable');
if isempty(common)
  error('apres:loadQuadpolSite:att', ...
    'no attenuator setting is shared by all four shots in %s', folder);
end
att = common(min(opt(opts, 'att_index', 1), numel(common)));
Q.attenuator = att;

for k = 1:numel(pols)
  o = opts;
  o.att_index = find(unique(B{k}.chirp_att, 'stable') == att, 1);
  P = apres.rangeProcess(B{k}, o);
  Q.(pols{k}) = P.spec;
  Q.time.(pols{k}) = P.timestamp;
  if k == 1
    Q.z = P.range; Q.fc = P.fc; Q.bandwidth = P.B; Q.ci = P.ci;
  elseif numel(P.range) ~= numel(Q.z) || P.fc ~= Q.fc
    error('apres:loadQuadpolSite:mismatch', ...
      '%s was recorded with different radar settings from HH', Q.files.(pols{k}));
  end
end

% SUB-SAMPLE CO-REGISTRATION. Each shot is a separate acquisition with
% the antennas moved, and the channels come out offset in range by a
% bin or two (0.2-3 bins on the GHOST sites). Misaligned speckle
% decorrelates, lowering |C_HHVV|. Measure each channel's offset against
% HH from the log-power envelope, then shift it by that amount.
if opt(opts, 'coregister', true)
  zi = Q.z > opt(opts, 'coreg_zmin', 40) & Q.z < opt(opts, 'coreg_zmax', 1000);
  for k = 2:numel(pols)
    lag = envelopeLag(Q.hh(zi), Q.(pols{k})(zi));
    Q.(pols{k}) = fracShift(Q.(pols{k}), -lag);
    Q.offset_bins.(pols{k}) = lag;
  end
end
[~, Q.site] = fileparts(folder);
end

% -------------------------------------------------------------------------
function lag = envelopeLag(a, b)
% Offset of b relative to a, in bins, from the peak of the cross-correlation
% of their log-power envelopes, refined to a fraction of a bin by fitting a
% parabola through the peak and its neighbours.
ea = log(movmean(abs(a).^2, 5)); ea = ea - mean(ea);
eb = log(movmean(abs(b).^2, 5)); eb = eb - mean(eb);
L = -10:10;
r = zeros(size(L));
for i = 1:numel(L)
  if L(i) >= 0, r(i) = sum(ea(1:end-L(i)) .* eb(1+L(i):end));
  else,         r(i) = sum(ea(1-L(i):end) .* eb(1:end+L(i))); end
end
[~, i] = max(r);
lag = L(i);
if i > 1 && i < numel(r)
  lag = lag + 0.5*(r(i-1) - r(i+1)) / (r(i-1) - 2*r(i) + r(i+1));
end
end

function y = fracShift(x, sh)
% Shift a complex profile by sh bins (fractional allowed) using the Fourier
% shift theorem: a delay is a linear phase ramp in the transform domain.
% The profile is oversampled (zero-padded FFT), so it is band-limited and
% this interpolation is exact up to edge effects.
x = x(:); N = numel(x);
k = [0:ceil(N/2)-1, -floor(N/2):-1]';
y = ifft(fft(x) .* exp(-2i*pi*k*sh/N));
end

function v = opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
