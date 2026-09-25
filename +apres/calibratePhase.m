function Q = calibratePhase(Q, opts)
%CALIBRATEPHASE Remove the antenna phase offsets between the four shots.
%
% Q = apres.calibratePhase(Q, opts)
%
% WHY. A quad-pol ApRES site is four separate shots, and the antennas are
% picked up and turned between them. Moving an antenna by a few cm changes
% the path length, and at a 1 m wavelength that is tens of degrees of
% phase. HH-VV phase at 20 m depth on the GHOST sites is 71, -135, -82 and
% 139 degrees - where fabric has had no depth to act, it should be ~0.
%
% WHY IT MUST BE FIXED ON THE CHANNELS. Every synthesized azimuth mixes
% all four channels, and the synthesized coherence obeys
% C(psi + 90) = conj(C(psi)) exactly (turning 90 degrees swaps H and V).
% So a constant phase cannot sit on the synthesized field - a fit cannot
% absorb it afterwards. It has to come off the channels first.
%
% THE MODEL. Each shot's phase is a transmit-antenna term plus a
% receive-antenna term. The field log (GHOST24_PpRES_Log_OZ.pdf) gives the
% layout: HH = Tx H / Rx H, HV = Tx H / Rx V, VH = Tx V / Rx H,
% VV = Tx V / Rx V. Relative to HH, with a = receive V-minus-H and
% b = transmit V-minus-H:
%     HV carries a,   VH carries b,   VV carries a + b
% Two measurements pin a and b down:
%   1. RECIPROCITY. The ice gives HV = VH, so the phase of <hv conj(vh)>
%      is a - b. HV and VH are nearly identical (|C| ~ 0.99), so this is
%      very precise.
%   2. ISOTROPIC FIRN. Near the surface the ice has no fabric to speak of,
%      so <hh conj(vv)> should have zero phase; its measured phase is
%      -(a + b).
%
% opts
%   .surf_m   ([10 40])   depth range treated as isotropic
%   .cross_m  ([40 800])  depth range used for the reciprocity phase
%
% Adds Q.phase_cal with fields a, b (rad), the two coherences behind them
% (so a weak estimate is visible downstream) and branch_note.

if nargin < 2, opts = struct(); end
sr = opt(opts, 'surf_m', [10 40]);
cr = opt(opts, 'cross_m', [40 800]);

s = Q.z >= sr(1) & Q.z <= sr(2);
c = Q.z >= cr(1) & Q.z <= cr(2);
x_co = sum(Q.hh(s) .* conj(Q.vv(s)));                   % ~ exp(-i(a+b))
x_cr = sum(Q.hv(c) .* conj(Q.vh(c)));                   % ~ exp(i(a-b))
apb = -angle(x_co);
amb =  angle(x_cr);
a = (apb + amb) / 2;
b = (apb - amb) / 2;
% A 180-DEGREE AMBIGUITY the data cannot break. (a + pi, b + pi) fits both
% measurements equally well; it flips the sign of HV and VH together,
% which MIRRORS the fabric axis about the H antenna (theta -> -theta).
% Only the antenna polarity - which way each dipole's R/L/U/D arrow points
% in the field log - decides it. Take the branch with the smaller offsets
% (antennas moved by centimetres, not half a wavelength), and record that
% a choice was made.
if abs(angle(exp(1i*a))) + abs(angle(exp(1i*b))) > ...
   abs(angle(exp(1i*(a+pi)))) + abs(angle(exp(1i*(b+pi))))
  a = a + pi; b = b + pi;
end
a = angle(exp(1i*a)); b = angle(exp(1i*b));

Q.hv = Q.hv * exp(-1i * a);
Q.vh = Q.vh * exp(-1i * b);
Q.vv = Q.vv * exp(-1i * (a + b));
Q.phase_cal = struct('a', a, 'b', b, ...
  'coh_surface', abs(x_co) / sqrt(sum(abs(Q.hh(s)).^2) * sum(abs(Q.vv(s)).^2)), ...
  'coh_reciprocity', abs(x_cr) / sqrt(sum(abs(Q.hv(c)).^2) * sum(abs(Q.vh(c)).^2)), ...
  'branch_note', 'smaller-offset branch chosen; the other mirrors theta about H');
end

function v = opt(o, f, d)
if isfield(o, f) && ~isempty(o.(f)), v = o.(f); else, v = d; end
end
