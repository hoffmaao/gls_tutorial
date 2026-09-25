%% Lesson 4 - Chi-square: do the data agree with our error model?
%
% Every error bar so far came from C_d, the noise we SAID the data has.
% On real data there is one survey and no Monte Carlo. Can C_d be checked
% from the data alone?
%
% Yes: if sample i has sigma_i, its residual should be about sigma_i, so
% the average of (r_i/sigma_i)^2 should be about 1. That average is the
% REDUCED CHI-SQUARE. fabricGLS and traveltimeFabricML report it as
% out.chi2_dof, and fabricGLS uses it to widen its error bars.

clear; close all; rng(4);

%% 1. The reduced chi-square
% Lesson 2's setup: phase vs depth with a sigma per sample.
g = 0.3001;
z = (10:10:1000)';
N = numel(z);
m_true = [0.4; 0.05];
G = [ones(N,1), g*z];
P = size(G, 2);                      % number of unknowns

sd_true = 0.1 + 0.4*(z/1000);        % the noise the data really has
d = G*m_true + sd_true .* randn(N,1);

w  = 1 ./ sd_true.^2;
m  = (G' * (w.*G)) \ (G' * (w.*d));
r  = d - G*m;
chi2 = sum(w .* r.^2)                % sum of (r_i / sigma_i)^2
dof  = N - P                         % degrees of freedom
chi2_dof = chi2 / dof                % should be near 1

% N - P, not N: fitting P unknowns absorbs some noise, so the residuals
% are slightly smaller than the noise. With a prior (lesson 5) the count
% becomes N minus the sum of the resolution diagonal, as in fabricGLS.

%% 2. Three scenarios
% (i)   the right sigma
% (ii)  sigma believed HALF its true size
% (iii) right sigma, wrong MODEL: the fabric strengthens with depth, so the
%       phase curves, and we still fit a line
d_curved = G*m_true + g*0.03*(z.^2/1000) + sd_true .* randn(N,1);

scen = {'(i) right sigma',  sd_true,   d;
        '(ii) sigma too small', sd_true/2, d;
        '(iii) model too simple', sd_true, d_curved};

figure('Name', 'Lesson 4a', 'Position', [100 100 1200 360]);
for s = 1:3
  sd_believed = scen{s,2};
  ds = scen{s,3};
  ws = 1 ./ sd_believed.^2;
  ms = (G' * (ws.*G)) \ (G' * (ws.*ds));
  rs = ds - G*ms;
  c2 = sum(ws .* rs.^2) / (N - P);
  fprintf('%-24s chi2/dof = %6.2f\n', scen{s,1}, c2);

  subplot(1,3,s);
  plot(z, rs ./ sd_believed, 'o', 'MarkerSize', 4); hold on;
  yline([-1 1], '--', 'Color', [0.4 0.4 0.4]); yline(0, 'Color', [0.4 0.4 0.4]);
  ylim([-8 8]); grid on;
  xlabel('depth [m]'); ylabel('r_i / \sigma_i');
  title(sprintf('%s: \\chi^2/dof = %.1f', scen{s,1}, c2));
end
% (ii): residuals too big everywhere. (iii): a pattern, a smile. chi2/dof
% cannot tell them apart; the residual plot can. Always look at both.

%% 3. What fabricGLS does with it: widen, never shrink
% If chi2/dof = 4 the residuals are twice what C_d claims, so multiply every
% error bar by sqrt(chi2/dof). From fabricGLS:
%
%     sig = sig * sqrt(max(chi2_dof, 1));
%
% max(..., 1): a better-than-expected fit does not earn tighter error bars.
% On a nonlinear problem it more often means a local minimum.
R = 3000;
for s = 2:3
  sd_believed = scen{s,2};
  ws = 1 ./ sd_believed.^2;
  C_M = inv(G' * (ws.*G));
  est = zeros(R,1); scaled = zeros(R,1);
  for k = 1:R
    if s == 2, dk = G*m_true + sd_true .* randn(N,1);
    else,      dk = G*m_true + g*0.03*(z.^2/1000) + sd_true .* randn(N,1); end
    mk = (G' * (ws.*G)) \ (G' * (ws.*dk));
    c2 = sum(ws .* (dk - G*mk).^2) / (N - P);
    est(k) = mk(2);
    scaled(k) = sqrt(C_M(2,2)) * sqrt(max(c2, 1));
  end
  fprintf('\n%s\n', scen{s,1});
  fprintf('  raw claimed sigma %.5f | scaled %.5f | measured scatter %.5f\n', ...
    sqrt(C_M(2,2)), median(scaled), std(est));
  fprintf('  mean estimate %.4f vs truth %.4f  (bias %.4f)\n', ...
    mean(est), m_true(2), mean(est) - m_true(2));
end

%% TAKEAWAYS
% - chi2/dof near 1: C_d and the model agree with the data. Far above 1:
%   C_d is too small or the model is missing something. Plot r./sigma.
% - Scaling fixes (ii). In (iii) it widens the error bar but cannot remove
%   the bias: no error bar fixes a wrong model. That is why fabricGLS lets
%   theta vary with depth instead of forcing it constant (lesson 7).
%
%% TRY THIS
% 1. Set sd_believed = 2*sd_true. What is chi2/dof? Should the error bars
%    shrink?
% 2. Add a column g*z.^2 to G in scenario (iii). What happens to chi2/dof
%    and the bias?
