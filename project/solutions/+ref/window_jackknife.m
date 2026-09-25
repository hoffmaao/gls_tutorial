function J = window_jackknife(p, W)
%WINDOW_JACKKNIFE Reference solution, step 5: error bars by resampling.
%
% J = ref.window_jackknife(p, W)
%
% WHY. The linearized error bars of step 4 FAIL the Monte Carlo test: all
% azimuths are synthesized from the same four channels, so their errors
% are correlated, and a diagonal C_d counts them as independent (lesson 3).
% Depth rows ~ one coherence window apart ARE nearly independent, so:
%
% THE JACKKNIFE. Drop one row, refit (starting from the full answer - no
% new grid search), repeat for every row. The spread of those n answers,
% scaled by (n-1)/n, estimates the variance of the full-data answer:
%     var_jk = (n-1)/n * sum_i (p_i - mean(p_i))^2
% Whatever correlation lives INSIDE a row is carried along automatically.
%
% theta is an axis: each replicate is folded onto the full answer's branch
% before differencing (never subtract raw axis angles).
%
% J fields: sigma [3 x 1] (theta0, delta0, ddelta), reps [3 x n]

rows = numel(W.z);
reps = zeros(3, rows);
for i = 1:rows
  k = true(rows, 1); k(i) = false;
  Wi = struct('psi', W.psi, 'z', W.z(k), 'zc', W.zc, 'C', W.C(k, :), 'sigma', W.sigma(k, :));
  pi_ = ref.fit_window(Wi, struct('p0', p));
  if pi_(3) * sign(p(3)) < 0 || abs(angle(exp(2i*(pi_(1) - p(1))))) > pi/2
    pi_ = [pi_(1) + pi/2; -pi_(2); -pi_(3)];     % undo the ddelta >= 0 flip
  end
  pi_(1) = p(1) + angle(exp(2i*(pi_(1) - p(1)))) / 2;
  pi_(2) = p(2) + angle(exp(1i*(pi_(2) - p(2))));
  reps(:, i) = pi_;
end
d = reps - mean(reps, 2);
J = struct('sigma', sqrt((rows - 1) / rows * sum(d.^2, 2)), 'reps', reps);
end
