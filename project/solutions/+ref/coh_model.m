function H = coh_model(psi, z, zc, p)
%COH_MODEL Reference solution, step 1: the single-column HH-VV coherence.
%
% H = ref.coh_model(psi, z, zc, p)
%
% psi [1 x Np] azimuths (rad), z [Nz x 1] depths (m), zc window centre (m),
% p = [theta0; delta0; ddelta]. Returns H [Nz x Np] with gamma = 1.
%
% mu = cos 2(psi - theta0), delta(z) = delta0 + ddelta*(z - zc), and
%   num = (1-mu^2)/2 + ((1+mu^2)/2) cos(delta) + i mu sin(delta)
%   den = 1 - ((1-mu^2)/2) (1 - cos(delta)),  floored at 0.05
% as in ptt.quadpolFabricLS.

mu  = cos(2 * (psi(:).' - p(1)));          % [1 x Np]
dl  = p(2) + p(3) * (z(:) - zc);           % [Nz x 1]
A   = (1 - mu.^2) / 2;
B   = (1 + mu.^2) / 2;
num = A + B .* cos(dl) + 1i * mu .* sin(dl);   % implicit expansion -> [Nz x Np]
den = max(1 - A .* (1 - cos(dl)), 0.05);
H   = num ./ den;
end
