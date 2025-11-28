%% From a Frequency-Domain Willems' Lemma to Data-Driven Predictive Control
% Example 2: Frequency-domain data-driven simulation using noisy data
% Date: October 30, 2025
% Author: T.J. Meijer - Eindhoven University of Technology
% Contact: t.j.meijer@tue.nl

% See the LICENSE file in the project root for full license information.
clc; 
clear all;
close all;
set(0,'defaultTextInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(0,'DefaultAxesFontSize', 14);
set(0,'DefaultTextFontSize', 14);

% Fix seed for reproducibility
rng(28102025);

% Add path to dependencies and functions
addpath('utils');
addpath('functions');

%% System description: Linearized four tank system
run('four_tank_system.m');

%% Parameters
M = 50; % Number of frequencies
L = 30; % Length of simulated future input-output trajectory
L0 = 3*nx; % Length of past input-output trajectory (Note: L0 >= lag(\Sigma))
u_ini = [diag([1,0.5])*1e2*ones(nu,1) zeros(nu,L0-1)]; % Initial input sequence
u = [-diag([1,2])*1e2*ones(nu,1) zeros(nu,L-1)]; % Future input sequence
E = nu; % Number of experiments
P = 120; % Number of periods to measure
rem_P = 20; % Remove first rem_Per periods
sigma_n = sqrt(0.1); % Noise covariance

%% Data generation (noisy)
% Excited frequencies
w = pi*[0:M-1]/M;

% Transfer function
G = @(w) C*((exp(1i*w)*eye(nx)-A)\B);

% Input excitation
ud = zeros(nu,2*M*P,E);
for e = 1:E
    ud(e,:,e) = 10*generate_multisine(w,P);
end

% Simulate to generate off-line data
yd = zeros(ny,2*M*P,E);
nd = zeros(ny,2*M*P,E);
for e = 1:E
    x = zeros(nx,1); % System state

    % Simulate system
    for i = 1:2*M*P
        nd(:,i,e) = sigma_n^2*randn(ny,1); % Measurement noise
        yd(:,i,e) = C*x + nd(:,i,e);
        x = A*x + B*ud(:,i,e);
    end
end

% Discard first rem_Per periods to get rid of the transient
ud = ud(:,2*M*rem_P+1:end,:);
yd = yd(:,2*M*rem_P+1:end,:);
nd = nd(:,2*M*rem_P+1:end,:);

P = P - rem_P; % Remaining number of periods

% Frequency-domain data
Ud = zeros(nu,M,E);
Yd = zeros(ny,M,E);
for e = 1:E
    Udhat = fft(ud(:,:,e),[],2)./(2*M*P);
    Ydhat = fft(yd(:,:,e),[],2)./(2*M*P);
    Ud(:,:,e) = Udhat(:,1:P:M*P);
    Yd(:,:,e) = Ydhat(:,1:P:M*P);
end

%% Persistence of excitation
% Compute frequency-domain input data matrix of depth L0+L+nx
Psi_u = Psi_matrix(L0+L+nx,Ud,w,E); 

% Rank test
if(rank(Psi_u) == nu*(L0+L+nx))
    display(['Input spectrum is PE of order L0+L+nx=',num2str(L0+L+nx),'.']);
else
    warning('Input spectrum is *not* PE of sufficient order.');
end

%% Generate data matrices for simulation
% Compute real-valued frequency-domain data matrices of depth L0+L
%   Using real-valued matrices to avoid having to impose G=(G0,G1,G1^*)
%   structure
Psi_u_real = Psi_matrix_real(L0+L,Ud,w,E);
Psi_y_real = Psi_matrix_real(L0+L,Yd,w,E);

%% Data-driven simulation
% Compute true input-output trajectory (model-based for verification)
[y_verif,~] = lsim(sys,[u_ini u],0:L0+L-1); 
y_verif = y_verif';

% Initial output sequence
y_ini = y_verif(:,1:L0);

% Data-driven simulation (Proposition 1)
[y,g] = dd_sim(Psi_u_real,Psi_y_real,u_ini,y_ini,u,L0,L);

% % Optional: Reconstruct G=(G0,G1,G1^*) and use it to compute y
% Psi_y = Psi_matrix(L0+L,Yd,w,E);
% T_Re = T_matrix(M,E);
% G = T_Re*g;
% y = Psi_y(ny*L0+1:end,:)*G;

%% Plot results
figure();
subplot(2,1,1);
plot([0:L-1],reshape(y,[2,L]),'LineWidth',2); hold on; grid on;
plot([-L0:L-1],y_verif,'--k','LineWidth',1.2); 
% legend('$y_1$','$y_2$','Model-based simulation','Location','northeast');
xlim([-L0 L-1])
xlabel('k'); ylabel('$y_k$');

%% Error computation
err_abs = norm(reshape(y,[2,L])-y_verif(:,L0+1:end));
err_rel = err_abs/norm(y_verif(:,L0+1:end));

display(['The absolute error is ',num2str(err_abs),'.']);
display(['The relative error is ',num2str(err_rel),'.']);

%% Data-driven simulation using time-domain data directly
% Compute Hankel matrix
H_u = Hankel_matrix(L0+L,ud);
H_y = Hankel_matrix(L0+L,yd);

% Data-driven simulation 
[y,g] = dd_sim(H_u,H_y,u_ini,y_ini,u,L0,L);

% Plot results
figure();
subplot(2,1,1);
plot([0:L-1],reshape(y,[2,L]),'LineWidth',2); hold on; grid on;
plot([-L0:L-1],y_verif,'--k','LineWidth',1.2); 
% legend('$y_1$','$y_2$','Model-based simulation','Location','northeast');
xlim([-L0 L-1])
xlabel('k'); ylabel('$y_k$');

% Error computation
err_abs = norm(reshape(y,[2,L])-y_verif(:,L0+1:end));
err_rel = err_abs/norm(y_verif(:,L0+1:end));

display(['The absolute error using time-domain data is ',num2str(err_abs),'.']);
display(['The relative error using time-domain data is ',num2str(err_rel),'.']);