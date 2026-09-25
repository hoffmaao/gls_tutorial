%% Lesson 2 - Weighted least squares: trust good samples more
%
% Birefringent ice (aligned crystals) splits a radar wave into two
% polarizations with slightly different speeds. Their phase difference
% grows with depth in a uniform column:
%
%     phi(z) = phi0 + g * dlam * z
%
% dlam is the FABRIC CONTRAST we want; g is a known constant. A straight
% line again, slope g*dlam.
%
% The catch: phi comes from the HH-VV COHERENCE, which fades with depth,
% and a low-coherence phase is noisy. Equal weights would let noisy deep
% samples pull as hard as clean shallow ones.
%
% phi is treated as unwrapped here to keep the fit linear. The real code
% never unwraps (lesson 6).

clear; close all; rng(2);

%% 1. The conversion constant g (as in ptt.traveltimeFabricML)
c       = 0.299792458e9;     % speed of light [m/s]
fc      = 750e6;             % radar centre frequency [Hz]
eps_bar = 3.171;             % average permittivity of ice
deps    = 0.034;             % permittivity anisotropy of a crystal
n_ice   = sqrt(eps_bar);
g = 2*pi*fc * 2 * (deps/(2*n_ice)) / c     % rad/m per unit dlam

%% 2. How noisy is a phase? The coherence says
% For coherence |C| from N independent looks, the phase standard deviation
% (Cramer-Rao bound, CRB) is
%
%     sd_phi = sqrt( (1 - |C|^2) / (2 * N * |C|^2) )
%
Cmag_examples = [0.95; 0.7; 0.4; 0.15];
N_look = 20;
sd_examples = sqrt((1 - Cmag_examples.^2) ./ (2*N_look*Cmag_examples.^2));
table(Cmag_examples, sd_examples, rad2deg(sd_examples), ...
  'VariableNames', {'coherence', 'sd_phi_rad', 'sd_phi_deg'})

% A few degrees for a clean sample, tens for a faded one. This formula
% sets the data errors in ptt.fabricGLS and the 'crb' weights in
% ptt.quadpolFabricLS.

%% 3. A synthetic column with fading coherence
z      = (10:10:1000)';
N      = numel(z);
dlam   = 0.05;
phi0   = 0.4;                   % phase at the surface [rad]
m_true = [phi0; dlam];

Cmag = 0.97 * exp(-z/450);                            % fades with depth
sd   = sqrt((1 - Cmag.^2) ./ (2*N_look*Cmag.^2));     % one sigma per sample

G = [ones(N,1), g*z];
d = G*m_true + sd .* randn(N,1);   % noise of each sample's own size

fprintf('phase sigma: %.2f rad at %d m, %.2f rad at %d m\n', sd(1), z(1), sd(end), z(end));

%% 4. Weights, one unknown first
% Measure each residual in units of its own sigma and minimise
%
%     sum( (r_i / sd_i)^2 ) = sum( w_i * r_i^2 ),   w_i = 1/sd_i^2
%
% The WEIGHT is 1 over the variance. Lesson 1's formulas with w inside:
%
%     slope = sum(w .* x .* d) / sum(w .* x.^2),  variance = 1 / sum(w .* x.^2)
%
w = 1 ./ sd.^2;
w(1:3)'                          % shallow: large
w(end-2:end)'                    % deep: tiny

%% 5. The weighted normal equations
% With W = diag(w):
%
%     (G' W G) m = G' W d,     C_M = inv(G' W G)
%
% No sigma^2 in front: the weights are the inverse variances.
W = diag(w);
m_wls = (G' * W * G) \ (G' * W * d)
C_wls = inv(G' * W * G);

% w .* G scales row i by w(i) without building W. Same result, less memory;
% the fabric code writes it this way.
m_wls_fast = (G' * (w .* G)) \ (G' * (w .* d))

%% 6. Two more routes to the same answer
% (a) WHITENING: divide each row of G and d by its sigma, so every residual
%     has sigma 1, then ordinary least squares. Every GLS in this tutorial
%     is "whiten, then OLS".
m_white = (G ./ sd) \ (d ./ sd)
% (b) lscov, MATLAB's built-in:
m_lscov = lscov(G, d, w)

%% 7. Ordinary least squares on the same data
% OLS estimates one average sigma from the residuals (as polyfit and fitlm do).
m_ols  = G \ d;
s2_hat = sum((d - G*m_ols).^2) / (N - 2);
C_ols  = s2_hat * inv(G'*G);                %#ok<MINV>

fprintf('\n             truth     OLS       WLS\n');
fprintf('dlam       %7.4f   %7.4f   %7.4f\n', dlam, m_ols(2), m_wls(2));

%% 8. Monte Carlo: which is better, and whose error bar is honest?
R = 3000;
dl_ols = zeros(R,1); dl_wls = zeros(R,1); s_ols = zeros(R,1);
for k = 1:R
  d_k = G*m_true + sd .* randn(N,1);
  mo = G \ d_k;
  dl_ols(k) = mo(2);
  Ck = (sum((d_k - G*mo).^2)/(N-2)) * inv(G'*G);          %#ok<MINV>
  s_ols(k) = sqrt(Ck(2,2));
  mw = (G' * (w.*G)) \ (G' * (w.*d_k));
  dl_wls(k) = mw(2);
end

fprintf('\n        measured scatter   claimed sigma   ratio\n');
fprintf('OLS     %10.5f        %10.5f       %.2f\n', std(dl_ols), median(s_ols), median(s_ols)/std(dl_ols));
fprintf('WLS     %10.5f        %10.5f       %.2f\n', std(dl_wls), sqrt(C_wls(2,2)), sqrt(C_wls(2,2))/std(dl_wls));
fprintf('WLS scatters %.1fx less than OLS, on the same data.\n', std(dl_ols)/std(dl_wls));

%% 9. Plot
figure('Name', 'Lesson 2', 'Position', [100 100 1200 380]);

subplot(1,3,1);
plot(Cmag, z, 'LineWidth', 1.5); set(gca, 'YDir', 'reverse');
xlabel('coherence |C|'); ylabel('depth [m]'); title('Coherence vs depth'); grid on;

subplot(1,3,2);
errorbar(d, z, sd, 'horizontal', 'o', 'MarkerSize', 3, 'CapSize', 0, ...
  'Color', [0.6 0.6 0.6], 'MarkerFaceColor', [0.3 0.3 0.3]); hold on;
plot(G*m_ols, z, 'LineWidth', 1.5);
plot(G*m_wls, z, '--', 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse'); xlabel('phase \phi [rad]'); ylabel('depth [m]');
title('Phase vs depth, \pm1\sigma'); grid on;
legend('data', 'OLS', 'WLS', 'Location', 'southwest');

subplot(1,3,3);
edges = linspace(min([dl_ols; dl_wls]), max([dl_ols; dl_wls]), 60);
histogram(dl_ols, edges); hold on; histogram(dl_wls, edges);
xline(dlam, 'k', 'LineWidth', 1.5);
xlabel('estimated \Delta\lambda'); ylabel('count');
title(sprintf('\\Delta\\lambda estimates, %d noise realizations', R)); legend('OLS', 'WLS', 'truth'); grid on;

%% TAKEAWAYS
% - Weight each sample by 1/variance. With the right variances this is the
%   best a linear fit can do.
% - OLS is unbiased but noisier, and its error bar is wrong: it assumed one
%   noise level.
% - This is why fabricGLS refuses phase data without |C|: a wrong phase
%   weight silently rotates the fabric axis.
%
%% TRY THIS
% 1. quadpolFabricLS caps its weights at 25. Add a cap here. Does the
%    estimate worsen? Why might real data want it anyway?
% 2. Set N_look = 5. Which numbers change, and by how much?
% 3. Compute sd with |C| 0.1 too high. Is the answer still centred on the
%    truth? Is the error bar still honest? (Lesson 4 catches this.)
