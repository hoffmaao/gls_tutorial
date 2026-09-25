%% Lesson 1 - Least squares
%
% A reflector at depth z returns after a two-way travel time t = s * z,
% with s = 2*n_ice/c the slowness. We measure t at many depths, each with
% random timing jitter. What is s, and how well is it known?
%
% Three ways to the same answer: try every s, then the formula, then the
% matrix form that every later lesson and the fabric code use.

clear; close all; rng(1);

%% 1. Synthetic data with a known answer
% Test an inversion on data whose truth you know before using it on data
% whose truth you don't.
c      = 0.299792458;        % speed of light [m per ns]
n_ice  = 1.78;               % refractive index of ice
s_true = 2*n_ice/c           % slowness [ns/m]

z     = (20:20:1000)';       % 50 depths [m]
N     = numel(z)
sigma = 3;                   % timing jitter [ns], the same for every sample

noise = sigma * randn(N, 1);
t     = s_true * z + noise;  % the measurements

[z(1:5), s_true*z(1:5), t(1:5)]      % depth, true time, measured time

%% 2. Try every s and keep the best
% For a trial slowness the residuals are measured minus predicted, and the
% misfit is their sum of squares (lesson 0).
s_grid = linspace(11.85, 11.90, 1001);
misfit = zeros(size(s_grid));
for k = 1:numel(s_grid)
  r = t - s_grid(k) * z;
  misfit(k) = sum(r.^2);
end
[~, k_best] = min(misfit);
s_grid_best = s_grid(k_best)

figure('Name', 'Lesson 1a', 'Position', [100 100 520 380]);
plot(s_grid, misfit, 'LineWidth', 1.5); hold on;
plot(s_grid_best, misfit(k_best), 'o', 'MarkerSize', 8, 'LineWidth', 1.5);
xline(s_true, '--', 'truth');
xlabel('trial slowness s [ns/m]'); ylabel('misfit  \Sigma r^2  [ns^2]');
title('Misfit vs trial slowness'); grid on;

% ptt.quadpolFabricLS finds the fabric axis this way: every angle on a
% grid, keep the smallest misfit. Slow, but immune to false minima
% (lesson 6).

%% 3. The formula
% The misfit is sum( (t_i - s z_i)^2 ). Its slope in s is
% -2 sum( z_i (t_i - s z_i) ); at the minimum the slope is zero, so
%
%     s = sum(z.*t) / sum(z.^2)
%
s_formula = sum(z .* t) / sum(z .^ 2)

%% 4. The error bar for one unknown
% s is a weighted sum of the t's, so the rules for adding random errors give
%
%     var(s) = sigma^2 / sum(z.^2)
%
% Deep samples (large z) make sum(z.^2) large and the error small.
sigma_s = sqrt(sigma^2 / sum(z.^2))

%% 5. Check the error bar by repeating the experiment
% An error-bar formula is a claim. Repeat the experiment with fresh noise
% and measure how much the answer scatters. This is Monte Carlo, and it is
% what fabric_anisotropy's test_fabric_gls does to fabricGLS.
R = 2000;
s_mc = zeros(R, 1);
for k = 1:R
  t_k = s_true * z + sigma * randn(N, 1);
  s_mc(k) = sum(z .* t_k) / sum(z .^ 2);
end
fprintf('claimed sigma %.6f   measured scatter %.6f   ratio %.2f\n', ...
  sigma_s, std(s_mc), sigma_s / std(s_mc));
% A ratio near 1 means the claimed error bar matches the scatter.

%% 6. Two unknowns
% Add an instrument delay t0:  t_i = 1 * t0 + z_i * s.
% A 2-D grid would still work; fabricGLS has 48 unknowns, so grids do not
% scale. Write the equations as a matrix instead (lesson 0):
t0_true = 25;
m_true  = [t0_true; s_true];
G = [ones(N,1), z];              % column 1 multiplies t0, column 2 multiplies s
d = G * m_true + sigma*randn(N, 1);   % d for data, from here on

G(1:4, :)                        % each row is one measurement's equation

%% 7. The normal equations
% G'*G and G'*d hold the sums from section 3:
G' * G          % [N, sum(z); sum(z), sum(z.^2)]
[N, sum(z); sum(z), sum(z.^2)]
G' * d          % [sum(d); sum(z.*d)]

% Section 3's  sum(z.^2) s = sum(z.*t)  becomes, with two unknowns,
%
%     (G'*G) * m = G'*d          the NORMAL EQUATIONS
%
m_normal = (G' * G) \ (G' * d)

% G\d solves the tall system directly and more accurately; polyfit does
% the same fit.
m_backslash = G \ d
p = polyfit(z, d, 1);            % returns [slope, intercept]
m_polyfit = [p(2); p(1)]

%% 8. The covariance matrix C_M
% Section 4's  sigma^2 / sum(z.^2)  becomes
%
%     C_M = sigma^2 * inv(G'*G)
%
% Its diagonal holds each unknown's variance (lesson 0, section 9).
C_M     = sigma^2 * inv(G' * G)      %#ok<MINV>
sigma_m = sqrt(diag(C_M))

% The off-diagonal: a steeper line must start lower, so the t0 and s
% errors anti-correlate.
rho = C_M(1,2) / (sigma_m(1) * sigma_m(2))

fprintf('t0 = %6.2f +- %.2f ns     (truth %.2f)\n', m_normal(1), sigma_m(1), t0_true);
fprintf('s  = %7.4f +- %.4f ns/m (truth %.4f)\n', m_normal(2), sigma_m(2), s_true);

%% 9. Check C_M by Monte Carlo
m_mc = zeros(2, R);
for k = 1:R
  d_k = G * m_true + sigma * randn(N, 1);
  m_mc(:, k) = G \ d_k;
end
measured = std(m_mc, 0, 2);
cc = corrcoef(m_mc(1,:), m_mc(2,:));
fprintf('\n            claimed sigma   measured scatter   ratio\n');
fprintf('t0 [ns]     %10.4f       %10.4f        %.2f\n', sigma_m(1), measured(1), sigma_m(1)/measured(1));
fprintf('s  [ns/m]   %10.6f       %10.6f        %.2f\n', sigma_m(2), measured(2), sigma_m(2)/measured(2));
fprintf('measured correlation %.2f vs claimed %.2f\n', cc(1,2), rho);

%% 10. Plot
figure('Name', 'Lesson 1b', 'Position', [100 100 1200 360]);

subplot(1,3,1);
plot(z, d, 'o', 'MarkerSize', 4); hold on;
plot(z, G*m_normal, 'LineWidth', 1.5);
xlabel('depth z [m]'); ylabel('two-way travel time [ns]');
title('Travel time vs depth'); legend('data', 'fit', 'Location', 'northwest');
grid on;

% Residuals should scatter evenly about zero with spread sigma. A pattern
% means the model is missing something (lesson 4).
subplot(1,3,2);
plot(z, d - G*m_normal, 'o', 'MarkerSize', 4); hold on;
yline([-sigma sigma], '--', 'Color', [0.4 0.4 0.4]);
yline(0, 'Color', [0.4 0.4 0.4]);
xlabel('depth z [m]'); ylabel('residual d - Gm [ns]');
title('Residuals'); grid on;

% Each Monte Carlo answer is a dot; the ellipse is where C_M puts 95% of
% them. Its tilt is the correlation.
subplot(1,3,3);
plot(m_mc(1,:), m_mc(2,:), '.', 'MarkerSize', 3, 'Color', [0.55 0.65 0.8]); hold on;
ang = linspace(0, 2*pi, 200);
[V, E] = eig(C_M);
ell = V * sqrt(E) * [cos(ang); sin(ang)] * sqrt(5.991);   % 95% for 2 unknowns
plot(m_true(1) + ell(1,:), m_true(2) + ell(2,:), 'LineWidth', 1.8);
xlabel('estimated t_0 [ns]'); ylabel('estimated s [ns/m]');
title('Estimates from 2000 noise realizations'); grid on;
legend('Monte Carlo', 'C_M, 95%', 'Location', 'northeast');

%% TAKEAWAYS
% - Least squares: the unknowns that minimise the sum of squared residuals.
%   The normal equations (G'G) m = G'd find them by setting the slope to zero.
% - C_M = sigma^2 inv(G'G) is the error bar, confirmed by Monte Carlo.
% - C_M assumed the same sigma for every sample, and independence between
%   samples. Lessons 2 and 3 drop those assumptions.
%
%% TRY THIS
% 1. Halve the depth range (z = (20:20:500)'). Which error bar grows more,
%    t0 or s? Why?
% 2. Use uniform noise: noise = sigma*sqrt(3)*(2*rand(N,1)-1). Is the
%    claimed sigma still right? (C_M only uses the variance.)
% 3. Let sigma grow with depth, sigma_i = 0.5 + 0.01*z_i, but keep one
%    sigma in C_M. Does the ratio stay near 1? That is lesson 2.
