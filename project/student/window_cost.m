function [cost, g, r] = window_cost(p, W)
%WINDOW_COST Step 2: the GLS misfit of one guess.
%
% [cost, g, r] = window_cost(p, W)
%
%   p     [theta0; delta0; ddelta]
%   W     one depth window: psi [1 x Np], z [Nz x 1] (each row's effective
%         depth), zc, C [Nz x Np] measured coherence, and
%         Wh [2Np x 2Np x Nz]: each row's whitening matrix
%   cost  chi-square misfit (a number)
%   g     best coherence scale for this p (real)
%   r     whitened residuals, one real column
%
% Each row's data are d_i = [Re C_i, Im C_i]' (2Np numbers). They are made
% from the same four channels, so their errors are correlated; Wh_i undoes
% that (lesson 3): Wh_i * d_i has independent, unit-variance entries.
% apres.coherenceField builds Wh from the exact covariance of the data.
%
% The model is g * H, with H from coh_model and h_i = [Re H_i, Im H_i]'.
% g enters linearly, so its best value is lesson 1's one-unknown formula
% in whitened units:
%     g = sum_i (Wh_i h_i)' (Wh_i d_i) / sum_i |Wh_i h_i|^2
% then g = max(g, 0): a negative scale would mean the model's pattern is
% upside down relative to the data, which no ice produces. The residuals
% and cost are
%     r = [Wh_1 (d_1 - g h_1); Wh_2 (d_2 - g h_2); ...],   cost = r' * r
%
% A tidy way to do every row at once: reshape d and h to [2Np x 1 x Nz]
% and use pagemtimes(W.Wh, d), which multiplies page i of Wh by page i of d:
%     d = reshape([real(W.C), imag(W.C)].', [], 1, Nz);
%
% Check:  check_step(2)

% TODO 2a: H from coh_model; d and h as [2Np x 1 x Nz].
H = [];

% TODO 2b: the whitened data and model, then g.
g = [];

% TODO 2c: r (one column) and cost.
r = [];
cost = [];

if isempty(H) || isempty(g) || isempty(r) || isempty(cost)
  error('project:todo', 'window_cost is not written yet - see the TODOs in project/student/window_cost.m');
end
end
