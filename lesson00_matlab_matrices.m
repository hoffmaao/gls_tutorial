%% Lesson 0 - The linear algebra you need, in MATLAB
%
% No background assumed. Run one section at a time (Ctrl+Enter, Cmd+Enter
% on a Mac) and read what MATLAB prints. Lines without a semicolon print
% their result on purpose.
%
% The whole tutorial in one sentence: write many equations as a table of
% numbers (a matrix) and let MATLAB find the numbers that satisfy them best.

clear; close all;

%% 1. A vector is a list of numbers
% A column vector runs top to bottom; ; starts a new row.
z = [100; 200; 300]          % three depths, in metres

% ' transposes: column to row and back.
z_row = z'

numel(z)                     % number of entries
size(z)                      % rows, columns
size(z_row)

%% 2. The dot operators act entry by entry
% A number times a vector scales every entry:
s = 11.9;                    % slowness, ns per metre
t = s * z

% .* .^ ./ act entry by entry. Plain * means something else (section 4).
a = [1; 2; 3];
b = [10; 20; 30];
a .* b                       % [1*10; 2*20; 3*30]
a .^ 2
a ./ b

%% 3. A sum of squares: the misfit
% The residual is measured minus predicted, one per sample:
t_measured  = [1195; 2375; 3574];
t_predicted = s * z;
r = t_measured - t_predicted

% The misfit squares each residual (so misses in both directions count)
% and adds them up:
misfit = sum(r.^2)

% r' * r is the same sum of squares: row times column multiplies matching
% entries and adds.
r' * r

%% 4. A matrix, and matrix times vector
G = [1 100;
     1 200;
     1 300]
size(G)                      % 3 rows, 2 columns

% With two unknowns, a delay t0 and slowness s, each depth gives
%     t_i = 1*t0 + z_i*s
% Stack the unknowns as m = [t0; s]. G * m does every equation at once:
% row i of G times m = G(i,1)*t0 + G(i,2)*s.
m = [25; 11.9];
G * m
1*25 + 100*11.9              % the first row, by hand

% Each row of G is one measurement; each column says how one unknown
% affects every measurement. G is the FORWARD OPERATOR, and every inverse
% problem here reads
%     d = G * m + noise        data = forward model of the unknowns + noise
% The job is to get m back from d.

%% 5. Building G from a vector
% ones(N,1) is a column of ones; [ , ] glues columns side by side.
z = (100:100:500)'
N = numel(z);
G = [ones(N,1), z]

%% 6. Solving equations: backslash
% Square system, exact solution:
%     2x + 1y = 5
%     1x + 3y = 10
A = [2 1; 1 3];
rhs = [5; 10];
xy = A \ rhs
A * xy                       % gives back [5; 10]

% With more equations than unknowns (the usual case) and noise, no m fits
% exactly. Backslash returns the m with the smallest sum of squared
% residuals: least squares, lesson 1.
d = [1210; 2400; 3590; 4780; 5960];     % noisy travel times at 5 depths
m_best = G \ d
residuals = d - G*m_best

%% 7. Identity, inverse, diag
I = eye(3)                   % ones on the diagonal: multiplying by it changes nothing
I * [4; 5; 6]

% inv(A) undoes A: the matrix version of 1/x. We use it for error bars,
% where the inverse itself is the answer. To solve, prefer A\b: same
% result, faster, more accurate.
inv(A) * A

% diag: vector -> diagonal matrix, or matrix -> its diagonal.
D = diag([1 4 9])
diag(D)

%% 8. Randomness: faking measurement noise
% randn draws from the bell curve, mean 0, standard deviation 1.
rng(0);                      % fix the random numbers so reruns match
noise = 3 * randn(100000, 1);
mean(noise)                  % ~0
std(noise)                   % ~3, the STANDARD DEVIATION
var(noise)                   % ~9, the VARIANCE = std^2

% Formulas are simpler in variances; results are quoted as standard
% deviations, which have the units of the measurement. "+- sigma" always
% means one standard deviation.

%% 9. Covariance: do two uncertain numbers wobble together?
e1 = randn(100000, 1);
e2 = randn(100000, 1);
x = e1;
y = -0.8*e1 + 0.6*e2;        % y partly copies MINUS e1

C = cov(x, y)                % 2x2 covariance matrix

% C(1,1), C(2,2): variances of x and y. C(1,2) = C(2,1): how they vary
% together. Divide by the two standard deviations for the CORRELATION,
% between -1 and 1:
rho = C(1,2) / sqrt(C(1,1) * C(2,2))   % ~-0.8: when x is high, y is low

figure('Name', 'Lesson 0', 'Position', [100 100 420 400]);
plot(x(1:3000), y(1:3000), '.', 'MarkerSize', 4);
axis equal; grid on; xlabel('x'); ylabel('y');
title(sprintf('Two correlated quantities, \\rho = %.2f', rho));

% Every inversion here returns its error bars as such a matrix, C_M: the
% square roots of its diagonal are the +- error bars, and the off-diagonal
% says which unknowns the data could not tell apart.

%% CHEAT SHEET
%   [1; 2; 3]        column vector            x'       transpose
%   .*  ./  .^       entry-by-entry           A * x    each row of A dotted with x
%   r' * r           sum of squares           A \ b    solve / least squares
%   eye(n)           identity                 inv(A)   undo A (use \ to solve)
%   diag(v)          vector -> diagonal       diag(A)  diagonal -> vector
%   randn            bell-curve noise         std/var  spread of a sample
%   cov(x, y)        covariance matrix        sqrt(diag(C))  error bars
