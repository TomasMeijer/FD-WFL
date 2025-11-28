%% From a Frequency-Domain Willems' Lemma to Data-Driven Predictive Control
% System description: Linearized four tank system
% Date: November 3, 2025
% Author: T.J. Meijer - Eindhoven University of Technology
% Contact: t.j.meijer@tue.nl

% See the LICENSE file in the project root for full license information.

%% Discrete-time system
% System matrices (source: J. Berberich, J. Köhler, M.A. Müller, and F. Allgöwer, 
% “Data-driven model predictive control with stability and robustness guarantees,”
% IEEE Trans. Autom. Control, vol. 66, no. 4, pp. 1702–1717, 2021.)
A = [0.921 0 0.041 0; 
    0 0.918 0 0.033; 
    0 0 0.924 0; 
    0 0 0 .937];
B = [0.017 0.001;
    0.001 0.023;
    0 0.061;
    0.072 0];
C = [1 0 0 0;
    0 1 0 0];

% System dimensions
[nx,nu] = size(B);
ny = size(C,1);

D = zeros(ny,nu);

sys = ss(A,B,C,[],-1);

% Transfer function
G = @(w) C*((exp(1i*w)*eye(nx)-A)\B);