function [y_fut,g] = dd_sim(Psi_u_real,Psi_y_real,u_ini,y_ini,u_fut,L0,L)
%DD_SIM Performs data-driven simulation
% 
% INPUTS:
%   Psi_u_real : Real-valued input data matrix (can be time-domain (Hankel) or frequency-domain)
%   Psi_y_real : Real-valued output data matrix (can be time-domain (Hankel) or frequency-domain)
%   u_ini : Initial input sequence
%   y_ini : Initial output sequence
%   u_fut : Future input sequence
%   L0 : Length of initial input-output trajectory
%   L : Length of future simulation
% 
% OUTPUTS:
%   y_fut : Simulated future output sequence
%   G : Corresponding latent variable g
% 
% Date: October 28, 2025
% Author: T.J. Meijer - Eindhoven University of Technology
% Contact: t.j.meijer@tue.nl

% See the LICENSE file in the project root for full license information.

%% Dimensions
nu = size(u_ini,1);
ny = size(y_ini,1);

%% Simulation
% Solve for G
g = pinv([Psi_u_real(1:nu*(L0+L),:); Psi_y_real(1:ny*L0,:)])*[u_ini(:); u_fut(:); y_ini(:)];

% Compute future output sequence
y_fut = Psi_y_real(ny*L0+1:end,:)*g;

end