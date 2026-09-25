function [cost, g, r] = window_cost(p, W)
%WINDOW_COST Step 2: the weighted misfit of one guess.
%
% [cost, g, r] = window_cost(p, W)
%
%   p     [theta0; delta0; ddelta]
%   W     one depth window: psi [1 x Np], z [Nz x 1], zc, C [Nz x Np]
%         measured coherence, sigma [Nz x Np] noise per real/imaginary part
%   cost  chi-square misfit (a number)
%   g     best coherence scale for this p (real, 0 to 1)
%   r     whitened residuals, one real column
%
% The model is g * H, with g the coherence loss. For fixed p, g enters
% linearly, so its best value is lesson 1's formula with lesson 2's weights:
%     w = 1 ./ sigma.^2
%     g = sum(w .* real(conj(H) .* C), 'all') / sum(w .* abs(H).^2, 'all')
% ptt.quadpolFabricLS does the same, so only the nonlinear unknowns are
% searched.
%
% g is real because C(psi + 90) = conj(C(psi)) for synthesized azimuths
% (turning the antennas 90 degrees swaps H and V). A constant phase cannot
% be absorbed here; it is removed from the channels first
% (apres.calibratePhase).
%
% Whitened residuals (lesson 2): divide each misfit by its sigma, then
% stack real and imaginary parts:
%     e = (C - g*H) ./ sigma;
%     r = [real(e(:)); imag(e(:))];
%     cost = r' * r
%
% Check:  check_step(2)

% TODO 2a: H from coh_model with g = 1.
H = [];

% TODO 2b: w and g.
g = [];

% TODO 2c: r and cost.
r = [];
cost = [];

if isempty(H) || isempty(g) || isempty(r) || isempty(cost)
  error('project:todo', 'window_cost is not written yet - see the TODOs in project/student/window_cost.m');
end
end
