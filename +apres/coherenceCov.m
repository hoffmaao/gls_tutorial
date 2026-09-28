function [Cd, C] = coherenceCov(M, N, psi)
%COHERENCECOV Covariance of the synthesized HH-VV coherences of one row.
%
% [Cd, C] = apres.coherenceCov(M, N, psi)
%
%   M    4x4 moments <s s'> of the channels (hh, vv, hv, vh), averaged
%        over N independent looks
%   psi  [1 x Np] azimuths, in the synthesis sense of apres.observables
%   C    [Np x 1] the coherences M gives at those azimuths
%   Cd   [2Np x 2Np] covariance of [Re C; Im C]
%
% Each coherence is C = a / sqrt(b c), with a = w_hh' M w_vv,
% b = w_hh' M w_hh, c = w_vv' M w_vv, and the w fixed weights for each
% azimuth. For complex Gaussian channels (Isserlis' theorem), bilinear
% forms alpha = B(p1,q1), beta = B(p2,q2) of the sample moments, with
% B(p,q) = p' M q, obey
%     E[d alpha conj(d beta)] = B(p1,p2) B(q2,q1) / N
%     E[d alpha d beta]       = B(p1,q2) B(p2,q1) / N
% Linearizing, dC = da/sqrt(bc) - (C/2)(db/b + dc/c). First order in 1/N:
% about 5% low at 10 looks (checked by simulation in lesson 9).

Np = numel(psi);
P1 = zeros(4, 3*Np); Q1 = P1;
for j = 1:Np
  c = cos(-psi(j)); s = sin(-psi(j)); cs = c*s;
  w_hh = [c*c; s*s; cs; cs];
  w_vv = [s*s; c*c; -cs; -cs];
  P1(:, j) = w_hh;      Q1(:, j) = w_vv;
  P1(:, Np+j) = w_hh;   Q1(:, Np+j) = w_hh;
  P1(:, 2*Np+j) = w_vv; Q1(:, 2*Np+j) = w_vv;
end
PP = P1.' * M * P1;  QQ = Q1.' * M * Q1;  PQ = P1.' * M * Q1;
K  = PP .* QQ.' / N;
Jp = PQ .* PQ.' / N;
v = diag(PQ);
a = v(1:Np); b = real(v(Np+1:2*Np)); cc = real(v(2*Np+1:end));
C = a ./ sqrt(b .* cc);
L = [diag(1 ./ sqrt(b .* cc)), diag(-C ./ (2*b)), diag(-C ./ (2*cc))];
G = L * K * L';
P = L * Jp * L.';
Cuu = real(G + P) / 2;  Cvv = real(G - P) / 2;  Cuv = imag(P - G) / 2;
Cd = [Cuu, Cuv; Cuv.', Cvv];
Cd = (Cd + Cd.') / 2;
end
