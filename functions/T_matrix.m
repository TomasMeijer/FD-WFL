function T = T_matrix(M,E)
%T_MATRIX Transformation matrix T_Re used to convert complex-valued
%conditions to real-valued by exploiting structure
% 
% INPUTS:
%   M : Number of frequencies (excluding negative frequencies)
%   E : Number of experiments/data sets
% 
% OUTPUTS:
%   T : Transformation matrix T_Re
% 
% Date: October 28, 2025
% Author: T.J. Meijer - Eindhoven University of Technology
% Contact: t.j.meijer@tue.nl

% See the LICENSE file in the project root for full license information.

T = (1/2)*blkdiag(2*eye(E),[eye((M-1)*E) -1i*eye((M-1)*E); eye((M-1)*E) 1i*eye((M-1)*E)]);

end

