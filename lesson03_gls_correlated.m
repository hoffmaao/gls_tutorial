%% Lesson 3 - Generalized least squares: correlated noise
%
% Lessons 1 and 2 assumed each sample's noise is independent of its
% neighbours'. Coherence is estimated by averaging over a depth window,
% and neighbouring windows share samples, so their errors move together.
%
% The DATA COVARIANCE MATRIX C_d records the variance of each sample and
% the covariance of each pair. Weighted least squares is the case where
% C_d is diagonal.

clear; close all; rng(3);

%% 1. Correlated noise
% Draw independent noise and smooth it with a running average of L
% samples, as a multilook window does.
N     = 200;                          % depth samples
L     = 10;                           % samples per averaging window
sigma = 0.15;                         % phase noise per sample [rad]

white  = randn(N + L - 1, 1);
smooth = conv(white, ones(L,1)/sqrt(L), 'valid');     % running sum / sqrt(L)
% sqrt(L), not L, keeps each sample's sigma at 1: only the correlation
% differs from independent noise.
e_corr  = sigma * smooth;
e_indep = sigma * randn(N, 1);
[std(e_indep), std(e_corr)]           % the same size

figure('Name', 'Lesson 3a', 'Position', [100 100 900 320]);
plot(e_indep, 'Color', [0.6 0.6 0.6]); hold on;
plot(e_corr, 'LineWidth', 1.4);
xlabel('depth sample'); ylabel('noise [rad]'); grid on;
legend('independent', sprintf('correlated over %d samples', L));
title('Noise realizations, \sigma = 0.15 rad');

%% 2. The data covariance matrix C_d
% Diagonal: sigma^2. Off-diagonal: samples k apart share L-k of their L
% inputs, so their covariance is sigma^2 (1 - k/L), zero beyond k = L.
% toeplitz builds a matrix that depends only on |i-j|.
lag  = (0:N-1)';
rho  = max(0, 1 - lag/L);
C_d  = sigma^2 * toeplitz(rho);

C_d(1:6, 1:6) / sigma^2               % top corner, as correlations

% Check against simulation: each column of E is one noise realization.
Rsim = 20000;
E = sigma * conv2(randn(N + L - 1, Rsim), ones(L,1)/sqrt(L), 'valid');
C_emp = (E * E') / Rsim;
fprintf('formula vs simulated, lag 0..3: [%s] vs [%s]\n', ...
  sprintf('%.4f ', C_d(1,1:4)), sprintf('%.4f ', C_emp(50,50:53)));

%% 3. The problem: phase vs depth (lesson 2's model)
g = 0.3001;                          % rad/m per unit dlam
z = (1:N)' * 2;                      % 2 m spacing
m_true = [0.4; 0.05];                % [phi0; dlam]
G = [ones(N,1), g*z];

%% 4. Generalized least squares
% Lesson 2 had (G' W G) m = G' W d with W = 1/variance. For a full
% covariance, 1/variance becomes the inverse matrix:
%
%     (G' inv(C_d) G) m = G' inv(C_d) d,       C_M = inv(G' inv(C_d) G)
%
% With a diagonal C_d this is lesson 2. lscov(G, d, C_d) does it built in.

%% 5. Whitening
% chol gives a lower-triangular Lc with Lc*Lc' = C_d: the recipe that turns
% independent unit noise into this correlated noise. Lc\ undoes it.
Lc = chol(C_d, 'lower');
e_white = Lc \ e_corr;
fprintf('\nafter whitening: std %.2f (want 1), lag-1 correlation %.2f (want 0)\n', ...
  std(e_white), sum(e_white(1:end-1).*e_white(2:end))/sum(e_white.^2));

% Whiten G and d the same way, then ordinary least squares:
d = G*m_true + e_corr;
Gw = Lc \ G;
dw = Lc \ d;
m_gls   = Gw \ dw
m_lscov = lscov(G, d, C_d)                    % agrees
C_gls   = inv(Gw' * Gw);                      %#ok<MINV> = inv(G' inv(C_d) G)

%% 6. Three treatments, compared by Monte Carlo
% (a) IGNORE the correlation: diagonal C_d only
% (b) FULL GLS with the true C_d
% (c) INFLATE: keep a diagonal C_d but multiply each variance by L, the
%     number of correlated samples. This is what fabricGLS does through
%     opts.n_looks / n_indep_z instead of building an N x N matrix.
C_a = sigma^2 * inv(G'*G);                    %#ok<MINV> (a)'s claimed error bar
C_c = L * C_a;                                % (c)'s claimed error bar

R = 3000;
dl = zeros(R, 2);                             % columns: (a)=(c), (b)
Ecol = sigma * conv2(randn(N + L - 1, R), ones(L,1)/sqrt(L), 'valid');
for k = 1:R
  d_k = G*m_true + Ecol(:, k);
  ma = G \ d_k;                               % (a) and (c): same estimate
  mb = (Lc \ G) \ (Lc \ d_k);                 % (b)
  dl(k, :) = [ma(2), mb(2)];
end

rows = {'(a) ignore correlation', '(b) full GLS', '(c) inflate variance'};
meas = [std(dl(:,1)); std(dl(:,2)); std(dl(:,1))];
clm  = sqrt([C_a(2,2); C_gls(2,2); C_c(2,2)]);
disp(table(meas, clm, clm./meas, 'RowNames', rows, ...
  'VariableNames', {'measured_scatter', 'claimed_sigma', 'ratio'}))

%% 7. Plot
figure('Name', 'Lesson 3b', 'Position', [100 100 1100 380]);
subplot(1,2,1);
imagesc(C_d(1:40, 1:40) / sigma^2); axis image; colorbar;
xlabel('sample j'); ylabel('sample i');
title('C_d / \sigma^2 (first 40 samples)');

subplot(1,2,2);
b = bar([meas, clm]); grid on;
set(gca, 'XTickLabel', {'(a) ignore', '(b) full GLS', '(c) inflate'});
ylabel('\sigma of \Delta\lambda');
legend('measured', 'claimed', 'Location', 'northoutside', 'Orientation', 'horizontal');
title('\sigma of \Delta\lambda by method');

%% TAKEAWAYS
% - Correlated noise barely changes the estimate but, ignored, makes the
%   claimed sigma about sqrt(L) too small.
% - Full GLS: whiten with chol, then least squares.
% - Inflating a diagonal C_d gets the size of the error bar right without
%   an N x N matrix; fabricGLS does this.
% - fabricGLS once assumed that 4-channel synthesized azimuths were ~4
%   independent numbers; its test measured that inflation made the error
%   bars 2.4x too large, so the option is off. Measure, don't assume.
%
%% TRY THIS
% 1. Set L = 1. What happens to the three rows of the table?
% 2. Set L = 30. Does the sqrt(L) rule for (a) still hold?
% 3. Inflate by 4 instead of L in (c). How wrong is the error bar?
