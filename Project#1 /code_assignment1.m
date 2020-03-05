% Assignment_1 Digital Filters
%
% author: Dimitra Ntzioni
% date: March 2017

clear all
close all

n = 500; % number of time steps

%% signal
A = sqrt(0.15)*randn(n,1); 
A = A - mean(A); 

for j=1:n
  x(j) = A(j)*sin(pi/8 *j + pi/6); % input to AR process
end

v = sqrt(0.32)*randn(n,1);
v = v - mean(v); % white noise

d = zeros(n,1); % initialize

% Desired signal
for i=1:n
  d(i) = x(i) + v(i);
end

%% channel
u = zeros(n,1); % initialize

u(1) = v(1);
u(2) = 0.25*u(1)+v(2);
for i=3:n
  u(i) = 0.25*u(i-1) - 0.12*u(i-2) + v(i);
end

x = reshape(x,[n,1]);
figure(1)
plot([d x u])
legend({'d(n)', 'x(n)', 'u(n)'})

%% FIR filter
a = [1 -0.25 0.12 ; -0.25 1.12 0 ; 0.12 -0.25 1];
b = [0.32 ; 0 ; 0];
r = a\b;
R = [r(1) r(2) r(3) ; r(2) r(1) r(2) ; r(3) r(2) r(1)]; % autocorrelation E[u u']
p = [0.32 ; 0 ; 0]; % cross correlation E[u d]

wo = R\p;
pT = reshape(p,[1,3]);
J = 0.395-pT/R*p;
l_max = max(eig(R));
mu_max = 2/l_max;
T = [u [0; u(1:n-1)] [0; 0; u(1:n-2)]]; 

y = T*wo; 

ey = d-y;

figure(2)
plot([d y])
legend({'d(n)', 'y(n)'})

figure(3)
plot([ey x])
legend({'ey(n)', 'x(n)'})
%% Steepest descent
w = [-1; -1; -1]; 
mu = 10;

wt = w;

for i=1:n
  w = w + mu*(p-R*w); % Adaptation steps
  wt=[wt w];
end

%% parameter error
figure(4)
we = (wt - wo*ones(1,n+1)).^2;
e = sqrt(sum(we));

semilogy(e);
xlabel('time step n');
ylabel('Parameter error');
title('Parameter error');

%% Denoising the data given

m= 1350;

load ('sound.mat');
load ('noise.mat');


U = [u [0; u(1:(length(u)-1))] [0; 0; u(1:(length(u)-2))]];
R_ = U'*U/length(u);
p_ = d'*U/length(u);
wo_ = p_/R_;
wo_ = reshape(wo_,[length(wo_),1]);
y_ = U*wo_;

ey_ = d - y_;

ey_ = reshape(ey_,[length(ey_),1]);
sound(ey_,fs);

l_max_ = max(eig(R_));
mu_max_ = 2/l_max_;

w_ = [-1; -1; -1]; 
mu_ = 0.1;

wt_ = w_;
pT_ = reshape(p_,[length(p_),1]);

for i=1:m
  w_ = w_ + mu_*(pT_-R_*w_); % Adaptation steps
  wt_=[wt_ w_];
end


figure(5)
we_ = (wt_ - wo_*ones(1,m+1)).^2;
e_ = sqrt(sum(we_));

semilogy(e_);
xlabel('time step n');
ylabel('Parameter error');
title('Parameter error');