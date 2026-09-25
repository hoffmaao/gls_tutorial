%% Lesson 5 - Priors and the resolution matrix
%
% The two polarizations arrive at slightly different times, and the delay
% accumulates with depth:
%
%     dtau(z) = k * integral of dlam from the surface to z
%
% We measure dtau at every depth and want dlam at every depth: as many
% unknowns as data. The data alone cannot determine them. A PRIOR, a
% stated belief about what profiles are reasonable, makes the problem
% solvable. This lesson rebuilds ptt.traveltimeFabricML and checks the
% result against it.

clear; close all; rng(5);

%% 1. The forward model: a running total as a matrix
cumsum([1; 2; 3; 4])

% The same as a lower triangle of ones: row i adds entries 1..i.
T = tril(ones(4))
T * [1; 2; 3; 4]

% With spacing dz the integral is dz times the running total, and k
% converts to nanoseconds: G = k * dz * T, and d = G*m + noise with m the
% whole dlam profile.
g  = 2*pi*750e6 * 2 * (0.034/(2*sqrt(3.171))) / 0.299792458e9;  % rad/m per dlam (lesson 2)
k  = g / (2*pi*750e6) * 1e9          % ns per metre per unit dlam
dz = 20;
z  = (dz:dz:1500)';
N  = numel(z)                        % 75 depths, 75 unknowns
G  = k * dz * tril(ones(N));

%% 2. A true profile and its data
% A band of stronger fabric around 700 m.
dlam_true = 0.03 + 0.05*exp(-((z - 700)/150).^2);
sigma_ns  = 0.05;                            % travel-time noise [ns]
d = G*dlam_true + sigma_ns*randn(N,1);

%% 3. The unregularized inversion
% G is square and invertible, so G\d solves exactly. But undoing a running
% total takes differences of neighbouring data, and differences of noisy
% numbers are mostly noise.
dlam_naive = G \ d;
fprintf('naive: rms error %.3f, vs the signal it is chasing (%.2f)\n', ...
  rms(dlam_naive - dlam_true), 0.05);

% This is an ILL-POSED problem: small changes in the data cause large
% changes in the answer. Inverting an integral is always like this.

%% 4. A prior, one number first
% A measurement says x = 5 +- 1; before measuring we believed x = 3 +- 2.
% Combine them by inverse-variance weights, as in lesson 2: the prior is
% one more measurement.
d1 = 5; s1 = 1;          % data
p1 = 3; sp = 2;          % prior
x_post = (d1/s1^2 + p1/sp^2) / (1/s1^2 + 1/sp^2)
s_post = sqrt(1 / (1/s1^2 + 1/sp^2))

%% 5. A prior for a profile
% The prior must also say how neighbouring depths relate. The fabric code
% uses (H_gp in fabricGLS and traveltimeFabricML)
%
%     C_m(i,j) = sp^2 * exp( -|z_i - z_j| / Lcorr )
%
% sp: how far dlam may stray from the prior mean. Lcorr: over what distance
% it may change.
m_prior = 0.03 * ones(N,1);
sp      = 0.05;
Lcorr   = 200;                               % metres
C_m     = sp^2 * exp(-abs(z - z') / Lcorr);  % z - z' is an N x N table of distances

% To understand a prior, draw profiles from it. chol turns independent
% noise into noise with covariance C_m (lesson 3).
Lm = chol(C_m + 1e-10*eye(N), 'lower');
draws = m_prior + Lm * randn(N, 5);
% The draws wander slowly at large scales and are rough at small ones:
% that is what the exponential covariance says.

%% 6. The posterior
% Minimise
%
%   (d - Gm)' inv(C_d) (d - Gm)  +  (m - m_prior)' inv(C_m) (m - m_prior)
%        misfit to the data            penalty for leaving the prior
%
% Setting the slope to zero (lesson 1) gives
%
%   A   = G' inv(C_d) G + inv(C_m)
%   m   = inv(A) * ( G' inv(C_d) d + inv(C_m) m_prior )
%   C_M = inv(A)
%
% the same shape as section 4. These lines are the core of
% traveltimeFabricML and of each Gauss-Newton step in fabricGLS.
Cdi = ones(N,1) / sigma_ns^2;           % inverse data variances
Cmi = inv(C_m);                         % inverse prior covariance
A   = G' * (Cdi .* G) + Cmi;
C_M = inv(A);
m_post  = C_M * (G' * (Cdi .* d) + Cmi * m_prior);
sig_post = sqrt(diag(C_M));
fprintf('with prior: rms error %.4f\n', rms(m_post - dlam_true));

%% 7. The resolution matrix
% A prior always yields an answer, even with no data. The resolution
% matrix
%
%     R = C_M * G' inv(C_d) G
%
% says how much of each estimate came from the data: estimate = R * truth
% + prior + noise. Row i a spike at i: depth i fully resolved. Row i zero:
% the prior decided. The diagonal is what fabricGLS and traveltimeFabricML
% return as out.resolution.
Rmat = C_M * (G' * (Cdi .* G));
res  = diag(Rmat);
% Each row sums to about 1 but is spread over several depths: the estimate
% is a local average of the truth. A small diagonal on a fine grid means
% "averaged over neighbours", not "unknown".
row = find(z == 700);
fprintf('resolution at 700 m: diagonal %.2f, row sum %.2f\n', res(row), sum(Rmat(row,:)));

%% 8. A data gap
% Remove the data between 900 and 1200 m (rows of G and d) and redo
% sections 6-7.
keep = ~(z > 900 & z < 1200);
Gk = G(keep, :); dk = d(keep); Cdk = Cdi(keep);
C_Mk  = inv(Gk' * (Cdk .* Gk) + Cmi);
m_k   = C_Mk * (Gk' * (Cdk .* dk) + Cmi * m_prior);
res_k = diag(C_Mk * (Gk' * (Cdk .* Gk)));

% traveltimeFabricML reports NaN where the resolution is below res_min,
% rather than an estimate that is mostly the prior:
%     m(res < res_min) = NaN;
res_min = 0.15;
m_report = m_k;
m_report(res_k < res_min) = NaN;
fprintf('with gap: %d of %d depths abstain\n', nnz(res_k < res_min), N);

%% 9. Check against ptt.traveltimeFabricML
% It integrates with the trapezoid rule (half of each end sample) rather
% than a running total. With that change the two must agree.
tutorial_setup();
Wt = zeros(N);
dzv = diff([0; z]);
for i = 2:N
  Wt(i,:) = Wt(i-1,:);
  Wt(i,i-1) = Wt(i,i-1) + dzv(i)/2;
  Wt(i,i)   = Wt(i,i)       + dzv(i)/2;
end
Gt = k * Wt;
C_Mt = inv(Gt' * (Cdi .* Gt) + Cmi);
m_mine = C_Mt * (Gt' * (Cdi .* d) + Cmi * m_prior);

out = ptt.traveltimeFabricML(d, z, struct('prior', [0.03 sp Lcorr], ...
  'sigma_ns', sigma_ns, 'res_min', 0));
fprintf('mine vs ptt.traveltimeFabricML: max difference %.2e\n', max(abs(m_mine - out.dlam)));

%% 10. Plots
figure('Name', 'Lesson 5a', 'Position', [100 100 1250 420]);

subplot(1,4,1);
plot(dlam_naive, z, 'Color', [0.6 0.6 0.6]); hold on;
plot(dlam_true, z, 'k', 'LineWidth', 2);
set(gca, 'YDir', 'reverse'); xlabel('\Delta\lambda'); ylabel('depth [m]');
title('Unregularized G\\d'); legend('G\\d', 'truth', 'Location', 'southeast'); grid on;

subplot(1,4,2);
plot(draws, z); hold on; plot(dlam_true, z, 'k', 'LineWidth', 2);
set(gca, 'YDir', 'reverse'); xlabel('\Delta\lambda'); title('5 draws from the prior'); grid on;

subplot(1,4,3);
fill([m_post - 2*sig_post; flipud(m_post + 2*sig_post)], [z; flipud(z)], ...
  [0.8 0.87 0.95], 'EdgeColor', 'none'); hold on;
plot(m_post, z, 'LineWidth', 1.5); plot(dlam_true, z, 'k', 'LineWidth', 2);
set(gca, 'YDir', 'reverse'); xlabel('\Delta\lambda'); title('Posterior \pm2\sigma');
legend('\pm2\sigma', 'posterior', 'truth', 'Location', 'southeast'); grid on;

subplot(1,4,4);
plot(Rmat(row,:), z, 'LineWidth', 1.5);
set(gca, 'YDir', 'reverse'); xlabel('weight'); grid on;
title('Row of R at 700 m');

figure('Name', 'Lesson 5b', 'Position', [100 100 800 420]);
subplot(1,2,1);
plot(res_k, z, 'LineWidth', 2.5, 'Color', [0.85 0.33 0.1]); hold on;
plot(res, z, '--', 'LineWidth', 1.5, 'Color', [0 0.45 0.74]);
xline(res_min, '--', 'res\_min');
set(gca, 'YDir', 'reverse'); xlabel('diag(R)'); ylabel('depth [m]');
title('Resolution'); legend('gap 900-1200 m', 'all data', 'Location', 'southeast'); grid on;

subplot(1,2,2);
plot(m_k, z, ':', 'LineWidth', 1.2); hold on;
plot(m_report, z, 'LineWidth', 2); plot(dlam_true, z, 'k', 'LineWidth', 1.5);
yregion(900, 1200, 'FaceColor', [0.9 0.9 0.9]);
set(gca, 'YDir', 'reverse'); xlabel('\Delta\lambda');
title('Posterior with the gap');
legend('raw posterior', 'reported', 'truth', 'Location', 'southeast'); grid on;

%% TAKEAWAYS
% - Inverting an integral amplifies noise. A prior (mean + covariance)
%   makes the problem solvable and states the assumption as numbers.
% - The posterior is lesson 2's weighted average with the prior as one
%   more measurement: A = G'inv(C_d)G + inv(C_m), C_M = inv(A).
% - The resolution matrix says how much of each answer came from the data.
%   Where it is low, report NaN; the fabric code does.
% - With a linear G these error bars are exact. Lesson 6 loses that.
%
%% TRY THIS
% 1. Set Lcorr to 50 and to 1000. Draw from the prior again. What changes
%    in the posterior and the resolution?
% 2. Set sp = 0.005. What happens to the bump at 700 m?
% 3. Move m_prior to 0.06. With the gap, which depths change most?
