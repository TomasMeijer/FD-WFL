function [H,F,G,h,f,g] = FreePC_prep(sys,T,Q,R,umax,umin,ymax,ymin)
%MPC_prep Constructs required matrices and vectors for MPC (compatible with quadprog)
% This implementation does not work for non-zero direct feedthrough matrix D    
%
% The cost function is given by 
%    (1/2)*x'*H*x 
% The equality constraints are given by
%    F*x = f
% The inequalitry constraints are given by
%    G*x <= g
% 
% INPUTS:
%   sys : State-space system
%   T : Prediction horizon
%   Q : Stage cost associated with states
%   R : Stage cost associated with inputs
%   umax : Upper bounds on u
%   umin : Lower bounds on u
%   ymax : Upper bounds on y
%   ymin : Lower bounds on y
%
% OUTPUTS:
%   H : Quadratic part of cost function
%   F : Equality constraint matrix
%   f : Equality constraint rhs vector
%   G : Inequality constraint matrix
%   g : Inequality constraint rhs vector
% 
% Date: November 12th, 2025
% Author: T.J. Meijer - Eindhoven University of Technology
% Contact: t.j.meijer@tue.nl

% See the LICENSE file in the project root for full license information.

% System matrices
A = sys.A;
B = sys.B;
C = sys.C;

% Dimensions
[nx,nu] = size(B);
ny = size(C,1);

% Construct H and h
H = 2*blkdiag(kron(eye(T),R),kron(eye(T),C'*Q*C));
h = @(uref,yref) [-2*kron(eye(T),R)*uref; -2*kron(eye(T),C'*Q)*yref];

% Construct F and f
F = zeros(T*nx,T*(nx+nu));
F(1:nx,1:nu) = -B;
F(1:nx,T*nu+(1:nx)) = eye(nx);
for i = 2:T
    F((i-1)*nx+1:i*nx,(i-1)*nu+(1:nu)) = -B;
    F((i-1)*nx+1:i*nx,T*nu+(i-2)*nx+(1:2*nx)) = [-A eye(nx)];
end
f = @(x0) [A*x0; zeros((T-1)*nx,1)];

% Construct G and g
G = [eye(nu*T) zeros(nu*T,nx*T);
    -eye(nu*T) zeros(nu*T,nx*T);
    zeros(ny*T,nu*T) kron(eye(T),C);
    zeros(ny*T,nu*T) -kron(eye(T),C)];
g = [kron(ones(T,1),umax); -kron(ones(T,1),umin); kron(ones(T,1),ymax); -kron(ones(T,1),ymin)];

end