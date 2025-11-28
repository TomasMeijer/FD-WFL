function [H,F,G,h,f,g] = FreePC_prep(Psi_u_real,Psi_y_real,T,Tbar,Q,R,lambda_sigma,lambda_G,umax,umin,ymax,ymin)
%FREEPC_prep Constructs required matrices and vectors for FreePC (compatible with quadprog)
%   
% The cost function is given by 
%    (1/2)*x'*H*x + h'*x
% The equality constraints are given by
%    F*x = f
% The inequalitry constraints are given by
%    G*x <= g
% 
% INPUTS:
%   Psi_u_real : Real-valued frequency-domain input data matrix
%   Psi_y_real : Real-valued frequency-domain output data matrix
%   T : Prediction horizon
%   Tbar : Length of past input-output trajectory
%   Q : Stage cost associated with states
%   R : Stage cost associated with inputs
%   lambda_sigma : Cost on 1-norm of sigma
%   lambda_G : Cost on 1-norm of G
%   umax : Upper bounds on u
%   umin : Lower bounds on u
%   ymax : Upper bounds on y
%   ymin : Lower bounds on y
%
% OUTPUTS:
%   H : Quadratic part of cost function
%   h : Linear part of cost function
%   F : Equality constraint matrix
%   f : Equality constraint rhs vector
%   G : Inequality constraint matrix
%   g : Inequality constraint rhs vector
% 
% Date: November 7th, 2025
% Author: T.J. Meijer - Eindhoven University of Technology
% Contact: t.j.meijer@tue.nl

% See the LICENSE file in the project root for full license information.

% Dimensions
nG = size(Psi_u_real,2);
nu = size(R,1);
ny = size(Q,1);

% Construct H and h (without regularization, but including sigma)
H = 2*blkdiag(kron(eye(T),R),kron(eye(T),Q),zeros(nG,nG),zeros(ny*Tbar,ny*Tbar)); % Standard cost on state, input, and G
h = @(uref,yref) [-2*kron(eye(T),R)*uref; -2*kron(eye(T),Q)*yref; zeros(nG,1); zeros(ny*Tbar,1)];

% Construct F (without regularization, but including sigma)
F = [zeros(nu*Tbar,nu*T) zeros(nu*Tbar,ny*T) Psi_u_real(1:nu*Tbar,:) zeros(nu*Tbar,ny*Tbar);
    -eye(nu*T) zeros(nu*T,ny*T) Psi_u_real(nu*Tbar+1:end,:) zeros(nu*T,ny*Tbar);
    zeros(ny*Tbar,nu*T) zeros(ny*Tbar,ny*T) Psi_y_real(1:ny*Tbar,:) -eye(ny*Tbar);
    zeros(ny*T,nu*T) -eye(ny*T) Psi_y_real(ny*Tbar+1:end,:) zeros(ny*T,ny*Tbar)]; % Standard prediction model
f = @(u_ini,y_ini) [u_ini; zeros(nu*T,1); y_ini; zeros(ny*T,1)];

% Construct G (without regularization, but including sigma)
G = [eye(nu*T) zeros(nu*T,ny*T) zeros(nu*T,nG) zeros(nu*T,ny*Tbar);
    -eye(nu*T) zeros(nu*T,ny*T) zeros(nu*T,nG) zeros(nu*T,ny*Tbar);
    zeros(ny*T,nu*T) eye(ny*T) zeros(ny*T,nG) zeros(ny*T,ny*Tbar);
    zeros(ny*T,nu*T) -eye(ny*T) zeros(ny*T,nG) zeros(ny*T,ny*Tbar);];
g = [kron(ones(T,1),umax); -kron(ones(T,1),umin); kron(ones(T,1),ymax); -kron(ones(T,1),ymin)];

% Add 1-norm regularization (by including slack variables t and s to implement
% 1-norm of sigma and G respectively)
H = blkdiag(H,zeros(ny*Tbar+nG));
h = @(uref,yref) [h(uref,yref); lambda_sigma*ones(ny*Tbar,1); lambda_G*ones(nG,1)];
F = [F,zeros(size(F,1),ny*Tbar+nG)];
G = [G,zeros(size(G,1),ny*Tbar+nG);
    zeros(ny*Tbar,(nu+ny)*T+nG) eye(ny*Tbar) -eye(ny*Tbar) zeros(ny*Tbar,nG); % sigma - t <= 0
    zeros(ny*Tbar,(nu+ny)*T+nG) -eye(ny*Tbar) -eye(ny*Tbar) zeros(ny*Tbar,nG); % -sigma - t <= 0
    zeros(nG,(nu+ny)*T) eye(nG) zeros(nG,2*ny*Tbar) -eye(nG); % G-s <= 0
    zeros(nG,(nu+ny)*T) -eye(nG) zeros(nG,2*ny*Tbar) -eye(nG); % -G-s <= 0
    zeros(ny*Tbar+nG,size(G,2)) -eye(ny*Tbar+nG)]; % t>=0 and s>= 0 
g = [g; zeros(3*(ny*Tbar+nG),1)];

end