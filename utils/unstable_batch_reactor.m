%% From a Frequency-Domain Willems' Lemma to Data-Driven Predictive Control
% System description: Unstable batch reactor with different outputs and
% pre-stabilizing controller
% Date: October 28, 2025
% Author: T.J. Meijer - Eindhoven University of Technology
% Contact: t.j.meijer@tue.nl

% See the LICENSE file in the project root for full license information.

%% Continuous-time system
% System matrices (source: G. C. Walsh and H. Ye, “Scheduling of networked
% control systems,” IEEE Control Syst. Mag., vol. 21, no. 1, pp. 57–65, 2001.)
Ac = [1.38 -0.2077 6.715 -5.676;
    -0.5814 -4.29 0 0.675;
    1.067 4.273 -6.654 5.893;
    0.048 4.273 1.343 -2.104];
Bc = [0 0;
    5.679 0;
    1.136 -3.146;
    1.136 0];
Cc = [1 0 1 -1; 0 1 0 0];

% System dimensions
[nx,nu] = size(Bc);
ny = size(Cc,1);

%% Discretization
Ts = 0.1; % Sampling time

A = expm(Ac*Ts);
B = (Ac\(A-eye(nx)))*Bc;
C = Cc;
D = zeros(ny,nu);
sys = ss(A,B,C,[],-1);

%% Transfer function
H = @(w) C*((exp(1i*w)*eye(nx)-A)\B);

%% Pre-stabilizing feedback controller
s = tf('s');
sys_K = [0 (2*s+2)/(s+0.05); (-5*s-8)/(s+0.05) 0]/1.84;
sys_K = ss(c2d(sys_K,Ts));
% z = tf('z');
% sys_K = [0 1.09*(z-0.9); -2.717*(z-0.840) 0]/(z-0.995);
% sys_K = ss(sys_K);
A_K = sys_K.A;
B_K = sys_K.B;
C_K = sys_K.C;
D_K = sys_K.D;

nd = 2;

%% Verify closed-loop stability
sys_cl = feedback(sys,sys_K);

% abs(eig(sys_cl.A))