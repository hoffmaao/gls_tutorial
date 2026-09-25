function Q = coregister(Q, opts)
%COREGISTER Align HV, VH and VV to HH in range, to a fraction of a bin.
%
% Q = apres.coregister(Q)
% Q = apres.coregister(Q, opts)
%
% Each shot is a separate acquisition with the antennas moved, and the
% channels come out offset in range by a bin or two (0.2-3 bins on the
% GHOST sites). Misaligned speckle decorrelates, lowering |C_HHVV|.
% Measure each channel's offset against HH from the log-power envelope,
% then shift it by that amount.
%
% opts: .coreg_zmin (40), .coreg_zmax (1000)  depths [m] the offset is
%       measured over
%
% Adds Q.offset_bins: measured (and removed) range offset of hv, vh, vv
% vs hh (bins).

if nargin < 2, opts = struct(); end
zi = Q.z > opt(opts, 'coreg_zmin', 40) & Q.z < opt(opts, 'coreg_zmax', 1000);
for p = {'hv', 'vh', 'vv'}
  lag = envelopeLag(Q.hh(zi), Q.(p{1})(zi));
  Q.(p{1}) = fracShift(Q.(p{1}), -lag);
  Q.offset_bins.(p{1}) = lag;
end
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
