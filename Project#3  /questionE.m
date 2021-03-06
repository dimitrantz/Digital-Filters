clear; close;

%% Initialize.
M = 150; % Number of filter coefficients.
delta = 100;
load music.mat
n = length(s);
s = s(:);
d = s;

u = zeros(size(d));
u(delta + 1:end) = d(1:end - delta);

%% Calculate correlations. Lecture 7 Page 3.
[r, lags] = xcorr(u, u, M + 1, 'unbiased');
% Autocorrelation vector (including r(0)).
r = r(lags >= 0);
% Autocorrelation matrix (needs r(0)).
R = toeplitz(r(1:end - 1));
% Autocorrelation vector (withour r(0).
rbar = r(2:end);
% Make sure that rbar is a column vector.
rbar = rbar(:);

%% Optimal forward linear prediction. Lecture 7 Page 4.
% Optimal Forward Predictor.
w_o = R \ rbar;
% Forward Prediction Error Power.
P_M = r(1) - rbar.' * w_o;

%% Call levinson function.
r = r(:);

% r is at least r(0), r(1), ..., r(M).
assert(M > 1);
assert(length(r) > M);

% Initialize.
P = [r(1), zeros(1, M)];
G = zeros(1, M);
L = zeros(M + 1, M + 1);
a = cell(1, M);
a{1} = 1;
L(1, 1) = a{1};
D = r(2);

% Main loop.
for m = 1:M
    G(m) = - D / P(m);
    a{m + 1} = [a{m}; 0] + G(m) * [0; a{m}(end:-1:1)];
    if (m < M)
        % We don't need to calculate D during the last iteration.
        D = a{m + 1}.' * r(m + 2:- 1:2);
    end
    P(m + 1) = P(m) * (1 - G(m) ^ 2);
    L(m + 1, 1:m + 1) = a{m + 1}(end:-1:1);
end


%% Cross correlation.
% Calculate the cross correlation between each signal b[i] and the desired
% signal s. Because b's size, n x (M+1), is too large (>4GB in current
% example) we caclulate it seperately for each coefficient and don't save
% the values afterwards.
fprintf('Calculating cross correlation vector p.\n');
p = zeros(M + 1, 1);
for i = 1:M + 1
   %
   fprintf('Coefficient %d/%d.\n', i, M + 1);

    b = filter(L(i, 1:i), 1, u);
    b(1:i) = u(1:i);

    [rb, lags] = xcorr(b, s, 1, 'unbiased');
    p(i) = rb(lags == 0);
end
fprintf('Done.\n');

%% Calculate the optimal parameters for the joint process estimator.
D = diag(P);
g_o = D \ p;

%%
fprintf('Calculating the predictor`s output.\n');
y = zeros(n, 1);
for i = 1:M + 1
   % fprintf('Coefficient %d/%d.\n', i, M + 1);

    b = filter(L(i, 1:i), 1, u);
    b(1:i) = u(1:i);
    y = y + g_o(i) * b;
end
fprintf('Done.\n');

e = s - y;

%% Play song.
fprintf('Playing song now. Press Enter to stop.\n');
sound(e, fs);
pause;
clear sound;