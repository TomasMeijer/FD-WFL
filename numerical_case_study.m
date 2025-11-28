%% From a Frequency-Domain Willems' Lemma to Data-Driven Predictive Control
% Numerical case study
% Date: October 31, 2025
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

% Add path to dependencies and functionsx_vb
addpath('utils');
addpath('functions');

% Set to 1 to avoid plotting results
no_plot = 0;

%% System description: Unstable batch reactor
run('unstable_batch_reactor.m');

%% Parameters
M = 40; % Number of frequencies
E = nu; % Number of experiments
P = 110; % Number of periods to measure
rem_P = 100; % Remove first rem_P periods
sigma_n = sqrt(0.3); % Noise covariance for off-line data collection

% FreePC/MPC parameters
T = 50;
Tbar = 3*nx;

% Simulation length
Lsim = 1e2;

%% Data generation (noisy)
% Excited frequencies
w = pi*[0:M-1]/M;

% Transfer function
G = @(w) C*((exp(1i*w)*eye(nx)-A)\B);

% Injected disturbance
dd = zeros(nu,2*M*P,E);
for e = 1:E
    dd(e,:,e) = 10*generate_multisine(w,P);
end

% Simulate to generate off-line data
ud = zeros(nu,size(dd,2),E);
yd = zeros(ny,size(dd,2),E);
for e = 1:E
    x = zeros(nx,1); % System state
    x_K = zeros(size(C_K,2),1); % Controller state

    % Simulate closed-loop system
    for i = 1:size(dd,2)
        n = sigma_n^2*randn(ny,1); % Measurement noise
        yd(:,i,e) = C*x + n;
        ud(:,i,e) = dd(:,i,e) + C_K*x_K - D_K*yd(:,i,e);
        x = A*x + B*ud(:,i,e);
        x_K = A_K*x_K - B_K*yd(:,i,e);
    end
end

% Discard first rem_Per periods to get rid of the transient
ud = ud(:,2*M*rem_P+1:end,:);
dd = dd(:,2*M*rem_P+1:end,:);
yd = yd(:,2*M*rem_P+1:end,:);

P = P - rem_P; % Remaining number of periods

% Averaging over periods
Syd_p = zeros(ny,nd,M,P);
Sud_p = zeros(nu,nd,M,P);
Hest_p = zeros(ny,nu,M,P);
for p = 1:P
    for m = 1:M
        for e = 1:E
            u_win = ud(:,(p-1)*2*M+1:p*2*M,e);
            d_win = dd(:,(p-1)*2*M+1:p*2*M,e);
            y_win = yd(:,(p-1)*2*M+1:p*2*M,e);
            U_win = fft(u_win,[],2)./(2*M);
            D_win = fft(d_win,[],2)./(2*M);
            Y_win = fft(y_win,[],2)./(2*M);

            Syd_p(:,:,m,p) = Syd_p(:,:,m,p) + Y_win(:,m)*D_win(:,m)';
            Sud_p(:,:,m,p) = Sud_p(:,:,m,p) + U_win(:,m)*D_win(:,m)';
        end
        Hest_p(:,:,m,p) = Syd_p(:,:,m,p)*inv(Sud_p(:,:,m,p));
    end
end
Hest = zeros(ny,nu,M);
Err = zeros(ny,nu,M);
Syd = zeros(ny,nd,M);
Sud = zeros(nu,nd,M);
for m = 1:M
    for p = 1:P
        Syd(:,:,m) = Syd(:,:,m) + Syd_p(:,:,m,p);
        Sud(:,:,m) = Sud(:,:,m) + Sud_p(:,:,m,p);
    end
    Hest(:,:,m) = Syd(:,:,m)*inv(Sud(:,:,m)); 
    % Hest(:,:,m) = H(w(m)); % Uncomment to use true FRF for debugging purposes
end

% Compute variance in estimates
Hvar = zeros(ny,nu,M);
for m = 1:M
    for p = 1:P
        for i = 1:ny
            for j = 1:nu
                Hvar(i,j,m) = Hvar(i,j,m) + (Hest_p(i,j,m,p)-Hest(i,j,m))'*(Hest_p(i,j,m,p)-Hest(i,j,m));
            end
        end
    end
end
Hvar = Hvar/(P*(P-1));

% Compute confidence intervals
conf_interval = 0.99;
Hmag_conf = zeros(ny,nu,M);
Hphase_conf = zeros(ny,nu,M);
for m = 1:M
    for i = 1:ny
        for j = 1:nu
            Hmag_conf(i,j,m) = sqrt(Hvar(i,j,m))*sqrt(2)*erfinv(conf_interval);
            Hphase_conf(i,j,m) = Hmag_conf(i,j,m)/abs(Hest(i,j,m));
        end
    end
end

% Frequency-domain data
Ud = zeros(nu,M,E);
Yd = zeros(ny,M,E);
for e = 1:E
    Ud(e,:,e) = 1;
    for m = 1:M
        Yd(:,m,e) = Hest(:,:,m)*Ud(:,m,e);
    end
end

% Plot FRF measurements
N = 1e3;
wplot = linspace(0,pi,N);
Htrue = zeros(ny,nu,N);
for i = 1:N 
    Htrue(:,:,i) = evalfr(sys,exp(1i*wplot(i)));
end

if ~no_plot
    % Bode plot without phase but with confidence intervals
    figure();
    h = tiledlayout(ny,nu);
    for i = 1:ny
        for j = 1:nu
            % Magnitude
            nexttile;
            plot(wplot,squeeze(20*log10(abs(Htrue(i,j,:)))),'k--','LineWidth',1);
            grid on; hold on;
            plot(w,squeeze(20*log10(abs(Hest(i,j,:)))),'b.');
            ylim_temp = ylim;
            fill([w,fliplr(w)],[squeeze(20*log10(abs(Hest(i,j,:))+Hmag_conf(i,j,:))).',fliplr(squeeze(20*log10(max(abs(Hest(i,j,:))-Hmag_conf(i,j,:),eps))).')],'b','facealpha',0.2,'linestyle','none');
            ylim([ylim_temp(1) ylim_temp(2)]);
            xline(pi,'k','LineWidth',1.1); 
            if j == 1
                ylabel(['Output ',num2str(i)],'FontSize',10);
            end
            if i == 1
                title(['Input ',num2str(j)],'FontSize',10);
            end
        end
    end
    xlabel(h,'Normalized frequency (rad/s)');
    ylabel(h,'Magnitude (db)');
end

if ~no_plot
    % Bode plot including phase
    figure();
    for i = 1:ny
        for j = 1:nu
            % Magnitude
            subplot(2*ny,nu,(i-1)*2*nu+j);
            plot(wplot,squeeze(20*log10(abs(Htrue(i,j,:)))),'k--','LineWidth',1);
            grid on; hold on;
            plot(w,squeeze(20*log10(abs(Hest(i,j,:)))),'b.');
            xline(pi,'k','LineWidth',1.1); 
            if j == 1
                ylabel('Magn. (db)');
            end
    
            % Phase
            subplot(2*ny,nu,(i-1)*2*nu+nu+j);
            plot(wplot,squeeze(rad2deg(unwrap(angle(Htrue(i,j,:))))),'k--','LineWidth',1);
            grid on; hold on;
            plot(w,squeeze(rad2deg(unwrap(angle(Hest(i,j,:))))),'b.');
            xline(pi,'k','LineWidth',1.1);
            if i == ny
                xlabel('Normalized frequency (rad/s)');
            end
            if j == 1
                ylabel('Phase (deg)');
            end
        end
    end
end

% Plot the error in FRF measurements
if ~no_plot
    figure();
    for i = 1:ny
        for j = 1:nu
            % Magnitude
            subplot(ny,nu,(i-1)*nu+j);
            plot(w,squeeze(20*log10(abs(Err(i,j,:)))),'r');
            grid on; hold on;
            xline(pi,'k','LineWidth',1.1); 
            if j == 1
                ylabel('Error (db)');
            end
        end
    end
end

%% Persistence of excitation
% Compute frequency-domain input data matrix of depth L0+L+nx
Psi_u = Psi_matrix(Tbar+T+nx,Ud,w,E); 

% Rank test
if(rank(Psi_u) == nu*(Tbar+T+nx))
    display(['Input spectrum is PE of order Tbar+T+nx=',num2str(Tbar+T+nx),'.']);
else
    warning('Input spectrum is *not* PE of sufficient order.');
end

%% Generate data matrices for FreePC
% Compute real-valued frequency-domain data matrices of depth L0+L
%   Using real-valued matrices to avoid having to impose G=(G0,G1,G1^*)
%   structure
Psi_u_real = Psi_matrix_real(Tbar+T,Ud,w,E);
Psi_y_real = Psi_matrix_real(Tbar+T,Yd,w,E);

%% Construct FreePC formulation 
% Reference
yref1 = [15; -10];
yref2 = [-5; 10];

% Compute steady-state input (model-based)
xrur = [A-eye(nx) B; C D]\[zeros(nx,1); yref1(1:ny)];
uref1 = xrur(nx+1:end);
xrur = [A-eye(nx) B; C D]\[zeros(nx,1); yref2(1:ny)];
uref2 = xrur(nx+1:end);

% Tuning parameters
Q = 1*eye(ny);
R = 1e-2*eye(nu);
lambda_sigma = 1e4;
lambda_G = 1e-2;

% Constraints
umax = [20; 20];
umin = [-10; -10];
ymax = [20; 20];
ymin = [-20; -20];

[Hfpc,Ffpc,Gfpc,hfpc,ffpc,gfpc] = FreePC_prep(Psi_u_real,Psi_y_real,T,Tbar,Q,R,lambda_sigma,lambda_G,umax,umin,ymax,ymin);

%% Frequency-domain data-driven predictive control (FreePC)
% Initialize initial input-output trajectory
y_ini0 = zeros(ny,Tbar);
u_ini0 = zeros(nu,Tbar);
u_ini = u_ini0;
y_ini = y_ini0;
x_ini = zeros(nx,1);

% Simulate closed-loop with FreePC
uk = zeros(nu,Lsim);
yk = zeros(ny,Lsim);
x0 = A*x_ini+B*u_ini(:,end);
xk = x0;
yref = kron(ones(T,1),yref1); 
uref = kron(ones(T,1),uref1);
cost = 0; % Store achieved cost
for k = 1:Lsim
    if k == Lsim/2
        yref = kron(ones(T,1),yref2);
        uref = kron(ones(T,1),uref2);
    end
    [xi,J,exitflag,output] = quadprog(Hfpc,hfpc(uref,yref),Gfpc,gfpc,Ffpc,ffpc(u_ini(:),y_ini(:)),[],[],[]);
    if exitflag ~= 1
        error('Quadprog failed');
    end
    uk(:,k) = xi(1:nu);
    
    % Simulate system response
    yk(:,k) = C*xk + D*uk(:,k); % Compute current measured output
    xk = A*xk + B*uk(:,k); % Compute next state

    % Update past input-output trajectory
    u_ini = [u_ini(:,2:end) uk(:,k)];
    y_ini = [y_ini(:,2:end) yk(:,k)];

    % Compute achieved cost
    cost = cost + (uk(:,k)-uref(1:nu))'*R*(uk(:,k)-uref(1:nu)) + (yk(:,k)-yref(1:ny))'*Q*(yk(:,k)-yref(1:ny));
end

if ~no_plot
    % Plot FreePC results
    figure();
    subplot(2,1,1);
    stairs([0:Lsim-1],uk','b','LineWidth',1.5); hold on;
    stairs([0;Lsim-1],[umin' umax'; umin' umax'],':','LineWidth',1.2);
    xlim([0,Lsim-1]); ylabel('$u_k$'); 
    ylim([-12,22]);
    subplot(2,1,2);
    stairs([0:Lsim-1],yk','b','LineWidth',1.5); hold on;
    stairs([0;Lsim-1],[ymin' ymax'; ymin' ymax'],':','LineWidth',1.2);
    stairs([0;Lsim/2;Lsim],[yref1';yref2';yref2'],'k','LineWidth',1);
    xlabel('$k$'); ylabel('$y_k$');
    xlim([0,Lsim-1]);
end

% Display achieved cost
display(['Achieved cost using FreePC: ',num2str(cost)]);

%% Model-based predictive control for comparison
% Construct MPC matrices
[Hmpc,Fmpc,Gmpc,hmpc,fmpc,gmpc] = MPC_prep(sys,T,Q,R,umax,umin,ymax,ymin);

% Simulate closed-loop with MPC
uk_mpc = zeros(nu,Lsim);
yk_mpc = zeros(ny,Lsim);
xk = x0;
yref = kron(ones(T,1),yref1); 
uref = kron(ones(T,1),uref1);
cost_mpc = 0; % Store achieved cost
for k = 1:Lsim
    if k == Lsim/2
        yref = kron(ones(T,1),yref2);
        uref = kron(ones(T,1),uref2);
    end
    [xi,J,exitflag,output] = quadprog(Hmpc,hmpc(uref,yref),Gmpc,gmpc,Fmpc,fmpc(xk),[],[],[]);
    if exitflag ~= 1
        error('Quadprog failed');
    end
    uk_mpc(:,k) = xi(1:nu);

    % Simulate system response
    yk_mpc(:,k) = C*xk + D*uk_mpc(:,k); % Compute current measured output
    xk = A*xk + B*uk_mpc(:,k); % Compute next state

    % Compute achieved cost
    cost_mpc = cost_mpc + (uk_mpc(:,k)-uref(1:nu))'*R*(uk_mpc(:,k)-uref(1:nu)) + (yk_mpc(:,k)-yref(1:ny))'*Q*(yk_mpc(:,k)-yref(1:ny));
end

% Plot MPC results (in FreePC plot)
subplot(2,1,1);
stairs([0:Lsim-1],uk_mpc','k--','LineWidth',1.2); 
subplot(2,1,2);
stairs([0:Lsim-1],yk_mpc','k--','LineWidth',1.2); 

% Plot MPC results (separately)
if ~no_plot
    figure();
    subplot(2,1,1);
    stairs([0:Lsim-1],uk_mpc','LineWidth',1.5); hold on;
    stairs([0;Lsim-1],[umin' umax'; umin' umax'],':','LineWidth',1.2);
    xlim([0,Lsim-1]); ylabel('$u_k$'); 
    ylim([-12,22]);
    subplot(2,1,2);
    stairs([0:Lsim-1],yk_mpc','LineWidth',1.5); hold on;
    stairs([0;Lsim-1],[ymin' ymax'; ymin' ymax'],':','LineWidth',1.2);
    stairs([0;Lsim/2;Lsim],[yref1';yref2';yref2'],'k','LineWidth',1);
    xlabel('$k$'); ylabel('$y_k$');
    xlim([0,Lsim-1]);
end

% Display achieved cost
display(['Achieved cost using MPC: ',num2str(cost_mpc)]);