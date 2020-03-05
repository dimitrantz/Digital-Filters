% clear/load
clear all
close all
load('speakerA.mat');
load('speakerB.mat');

%% Init

n = length(u);
coefficients = 6600; 
 

%% Wiener
fprintf('Wiener \n');
% init
P = xcorr(d, u);
P = P/length(u);
P = P(n:n+coefficients-1);
r = var(u)*autocorr(u, coefficients - 1); 
R = toeplitz(r);
wo = R \ P; 

% calculations
y = zeros(n, 1);
for i = coefficients + 1:n
     y(i) = wo' * u(i:-1:(i - coefficients + 1));
end
e = d - y; 
fprintf('sound: \n');
sound(e, fs);

%% LMS
fprintf('LMS \n');

% init
mu = 0.0009;
w = zeros(coefficients, 1);
e = zeros(n, 1);
y = zeros(n, 1);

% calculations
for i = (coefficients + 1):n
    y(i) = w' * u(i: - 1:(i - coefficients + 1));
    e(i) = d(i) - y(i);
    w = w + mu * e(i) * u(i:-1:(i - coefficients + 1));
end
fprintf('sound: \n');
sound(e, fs);

%% normalized LMS
fprintf('normalizes LMS \n');
% init
mu = 0.5; 
a = 100;
w = zeros(coefficients, 1);
e = zeros(n, 1);
y = zeros(n, 1);

for i = (coefficients + 1):n
    y(i) = w'*u(i: - 1:(i - coefficients + 1));
    e(i) = d(i) - y(i);
    w = w + mu * e(i) * u(i:-1:(i - coefficients + 1)) / (a + u(i:-1:(i - coefficients + 1))' * u(i:-1:(i - coefficients + 1)));
end
fprintf('sound: \n');
sound(e, fs);

%% RLS
fprintf('RLS \n');
% init
coefficients = 1000;
l = 1;
de = 0.005;
P = (1 / de) * eye(coefficients, coefficients);
w = zeros(coefficients, 1);
e = zeros(n, 1);
y = zeros(n, 1);


% calculations
for i = (coefficients + 1):n
    y(i) = w' * u(i:-1:(i - coefficients + 1));
    k = ((l^-1) * P * u(i:-1:i - coefficients + 1) / (1 + (l^-1) * u(i:-1:i - coefficients + 1)' * P * u(i:-1:(i - coefficients + 1))));
    e(i) = d(i) - y(i);
    w = w + k * e(i);
    P = (l^-1) *P - (l^-1) * k * u(i:-1:(i - coefficients + 1))' * P;
end
fprintf('sound: \n');
sound(e, fs)