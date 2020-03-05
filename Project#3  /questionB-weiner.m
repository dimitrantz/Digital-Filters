clear; close;

%% Initialize.
% Create the white noise.
sigmav2 = 0.54;  % Noise variance.
n = 250000;
v = sqrt(sigmav2) * randn(n, 1);
v = v - mean(v);

% Create the periodical interference.
f_o = 1 / 4;
phi = pi / 2;
A = 4.2;
i = 1:n;
x(i, 1) = A * (sin(2 * pi * f_o * i + phi) + cos(4 * pi * f_o * i + phi) + cos(7 * pi * i + phi / 3));

s = x + v;
delta = 10;  % delay.
M = 100;  % Number of filter coefficients.

% u(n) = s(n - delta).
u = zeros(size(s));
u(delta + 1:end) = s(1:end - delta);

%% Calculate correlations. Lecture 7 Page 3.
[r, lags] = xcorr(u, u, M + 1, 'unbiased');
r = r(lags >= 0);  % Autocorrelation vector (including r(0)).
R = toeplitz(r(1:end - 1));  % Autocorrelation matrix (needs r(0)).
rbar = r(2:end);  % Autocorrelation vector (withour r(0).
rbar = rbar(:);  % Make sure that rbar is a column vector.

%% Wiener filter coefficients.
w_o = R \ rbar;
