function J = window_jackknife(p, W)
%WINDOW_JACKKNIFE Step 5: error bars by resampling.
%
% J = window_jackknife(p, W)
%
% The step-4 error bars cannot be trusted (guide, Monte Carlo test): the
% 18 azimuths are synthesized from 4 channels, so their errors are
% correlated (lesson 3), and a sigma per entry treats them as independent.
% chi2/dof far below 1 says that sigma is wrong as well.
%
% Depth rows one coherence window apart are nearly independent. Leave one
% row out, refit from p, and repeat for each row. With n rows and
% replicate answers p_1..p_n:
%     sigma_jk = sqrt( (n-1)/n * sum_i (p_i - mean(p))^2 )
%
% For each row i:
%   1. Wi = W without row i (drop it from z, C and sigma)
%   2. q = fit_window(Wi, struct('p0', p))
%   3. theta0 is an axis. If theta0 moved by more than 45 degrees from p,
%      the refit landed on the twin (ddelta >= 0 is enforced), so switch
%      to q = [q(1) + pi/2; -q(2); -q(3)]. Then fold onto p's branch:
%          q(1) = p(1) + angle(exp(2i*(q(1) - p(1)))) / 2
%          q(2) = p(2) + angle(exp(1i*(q(2) - p(2))))
%
% J: struct with sigma [3 x 1], reps [3 x n]
%
% Check:  check_step(5)

% TODO: the leave-one-out loop, then sigma.
error('project:todo', 'window_jackknife is not written yet - see the TODOs in project/student/window_jackknife.m');
end
