%% From a Frequency-Domain Willems' Lemma to Data-Driven Predictive Control
% Example 3: Frequency-domain data-driven LQR
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
C = eye(nx); % Assume full state information
sys = ss(sys.A,sys.B,C,[],-1);
ny = nx;

%% Parameters
M = 50; % Number of frequencies
E = nu; % Number of experiments

% LQR cost matrices
Qlqr = eye(nx);
Rlqr = eye(nu);

%% Data generation (noisy)
% Excited frequencies
w = pi*[0:M-1]/M;

% Transfer function
G = @(w) C*((exp(1i*w)*eye(nx)-A)\B);

% Input spectrum
Ud = zeros(nu,M,E);
for e = 1:E
    Ud(e,:,e) = 1;
end

% Output spectrum
Yd = zeros(ny,M,E);
for m = 1:M
    for e = 1:E
        Yd(:,m,e) = G(w(m))*Ud(:,m,e);
    end
end

%% Persistence of excitation
% Compute frequency-domain input data matrix of depth L0+L+nx
Psi_u = Psi_matrix(nx+1,Ud,w,E); 

% Rank test
if(rank(Psi_u) == nu*(nx+1))
    display(['Input spectrum is PE of order nx+1=',num2str(nx+1),'.']);
else
    warning('Input spectrum is *not* PE of sufficient order.');
end

%% Data-driven LQR
% Generate data matrices for LQR
Psi_u = Psi_matrix(1,Ud,w,E);
Psi_x = Psi_matrix(2,Yd,w,E);

% Select X0, X1, U as in Proposition 3
X0 = Psi_x(1:nx,:);
X1 = Psi_x(nx+1:end,:);
U = Psi_u;

% Solve DARE using data-driven LMIs (Proposition 3)
P = sdpvar(nx);
lmi = [X0'*P*X0 - X1'*P*X1 - X0'*Qlqr*X0-U'*Rlqr*U <= -1e-7,P>=1e-7];
diagn = optimize(lmi,-trace(P),sdpsettings('solver','mosek'));
P = value(P);

% Construct LQR gain
N = null(X0'*P*X0 - X1'*P*X1 - X0'*Qlqr*X0-U'*Rlqr*U);
Xm_pinv = N*pinv(X0*N);
K = real(U*Xm_pinv);

%% Model-based LQR for reference
Pbar = dare(A,B,Qlqr,Rlqr);
Kbar = -(Rlqr+B'*Pbar*B)\B'*Pbar*A;

%% Error computation
err_abs_P = norm(P-Pbar);
err_rel_P = err_abs_P/norm(Pbar);
err_abs_K = norm(K-Kbar);
err_rel_K = err_abs_K/norm(Kbar);

display(['The absolute error \|P-Pbar\| is ',num2str(err_abs_P),'.']);
display(['The relative error \|P-Pbar\|/\|Pbar\| is ',num2str(err_rel_P),'.']);

display(['The absolute error \|K-Kbar\| is ',num2str(err_abs_K),'.']);
display(['The relative error \|K-Kbar\|/\|Kbar\| is ',num2str(err_rel_K),'.']);