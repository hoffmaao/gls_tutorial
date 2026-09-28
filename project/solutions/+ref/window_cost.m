function [cost, g, r] = window_cost(p, W)
%WINDOW_COST Reference solution, step 2: GLS misfit with g solved exactly.
%
% [cost, g, r] = ref.window_cost(p, W)
%
% W.psi [1 x Np], W.z [Nz x 1], W.zc, W.C [Nz x Np] complex data,
% W.Wh [2Np x 2Np x Nz] whitening matrix of each row (apres.coherenceField).
%
% Each row's data vector is d_i = [Re C_i, Im C_i]', and Wh_i * d_i has
% independent unit-variance entries (lesson 3). The model is g * H(p):
%     g = sum_i (Wh_i h_i)'(Wh_i d_i) / sum_i |Wh_i h_i|^2
% (lesson 1's one-unknown formula in whitened units), clamped at 0: a
% negative g would mean the pattern is upside down, which no ice makes.
%     r = stacked Wh_i (d_i - g h_i),   cost = r' * r.

H = ref.coh_model(W.psi, W.z, W.zc, p);
Nz = numel(W.z);
d = reshape([real(W.C), imag(W.C)].', [], 1, Nz);   % [2Np x 1 x Nz]
h = reshape([real(H),   imag(H)].',   [], 1, Nz);
dw = pagemtimes(W.Wh, d);
hw = pagemtimes(W.Wh, h);
g = sum(hw .* dw, 'all') / sum(hw.^2, 'all');
g = max(g, 0);                 % a coherence scale cannot be negative
r = reshape(dw - g * hw, [], 1);
cost = r' * r;
end
