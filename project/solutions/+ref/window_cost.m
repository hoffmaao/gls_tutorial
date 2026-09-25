function [cost, g, r] = window_cost(p, W)
%WINDOW_COST Reference solution, step 2: weighted misfit with gamma solved exactly.
%
% [cost, g, r] = ref.window_cost(p, W)
%
% W.psi [1 x Np], W.z [Nz x 1], W.zc, W.C [Nz x Np] complex data,
% W.sigma [Nz x Np] noise per real/imaginary part.
%
% The model is g * H(p) with g a real scale (the coherence loss). g enters
% linearly, so for fixed p the weighted least-squares g is a one-unknown
% problem (lesson 1, section 3, with weights):
%     g = sum(w .* real(conj(H) .* C)) / sum(w .* |H|^2),   w = 1 ./ sigma.^2
% g is real because the synthesized field obeys C(psi+90) = conj(C(psi)):
% a complex g would come out real anyway (see apres.calibratePhase).
% r is the WHITENED residual (lesson 2): real and imaginary parts of
% (C - g H) ./ sigma stacked into one real column. cost = r' * r.

H = ref.coh_model(W.psi, W.z, W.zc, p);
w = 1 ./ W.sigma.^2;
g = sum(w .* real(conj(H) .* W.C), 'all') / sum(w .* abs(H).^2, 'all');
e = (W.C - g * H) ./ W.sigma;
r = [real(e(:)); imag(e(:))];
cost = r' * r;
end
