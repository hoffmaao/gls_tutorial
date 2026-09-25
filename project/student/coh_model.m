function H = coh_model(psi, z, zc, p)
%COH_MODEL Step 1: the forward model.
%
% H = coh_model(psi, z, zc, p)
%
%   psi  [1 x Np] antenna azimuths [rad], a row
%   z    [Nz x 1] depths in the window [m], a column
%   zc   window centre depth [m]
%   p    [theta0; delta0; ddelta]: fabric axis azimuth [rad], phase
%        difference at zc [rad], and its rate with depth [rad/m]
%   H    [Nz x Np] complex coherence, one row per depth, one column per azimuth
%
% The model (lesson 6; ptt.quadpolFabricLS):
%   mu(psi)  = cos(2*(psi - theta0))
%   delta(z) = delta0 + ddelta*(z - zc)
%   A = (1 - mu.^2)/2,   B = (1 + mu.^2)/2
%   num = A + B.*cos(delta) + 1i*mu.*sin(delta)
%   den = max(1 - A.*(1 - cos(delta)), 0.05)
%   H   = num ./ den
%
% mu is a row (per azimuth) and delta a column (per depth). MATLAB expands
% a row combined with a column into a full table, so num and den come out
% [Nz x Np] without loops. Use psi(:).' and z(:) to fix the shapes.
%
% Check:  check_step(1)

% TODO 1a: mu, a row.
mu = [];

% TODO 1b: delta, a column.
dl = [];

% TODO 1c: A, B, num, den, H.
H = [];

if isempty(mu) || isempty(dl) || isempty(H)
  error('project:todo', 'coh_model is not written yet - see the TODOs in project/student/coh_model.m');
end
end
