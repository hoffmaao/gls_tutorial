function E = window_errors(p, W)
%WINDOW_ERRORS Step 4: linearized error bars and the chi-square check.
%
% E = window_errors(p, W)
%
% 1. J: the Jacobian at the answer, as in step 3 (same nudges h).
% 2. The residuals are whitened, so lesson 1's C_M = sigma^2 inv(G'G) is
%        C_M = inv(J' * J)
% 3. Chi-square (lesson 4): each real entry of r is one datum; the fit used
%    3 parameters plus g:
%        dof = numel(r) - 4
%        chi2_dof = cost / dof
% 4. Widen, never shrink (lesson 4; ptt.fabricGLS):
%        sigma = sqrt(diag(C_M)) * sqrt(max(chi2_dof, 1))
%
% E: struct with sigma [3 x 1], C_M [3 x 3], chi2_dof, dof, g
%
% Check:  check_step(4), then the Monte Carlo test in the guide.

% TODO 4a: cost, g, r at p; the Jacobian J.
% TODO 4b: C_M, dof, chi2_dof, sigma.
error('project:todo', 'window_errors is not written yet - see the TODOs in project/student/window_errors.m');

E = struct('sigma', sigma, 'C_M', C_M, 'chi2_dof', chi2_dof, 'dof', dof, 'g', g); %#ok<UNRCH>
end
