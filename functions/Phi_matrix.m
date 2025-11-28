function Phi = Phi_matrix(L,V,w,E)
%PSI_MATRIX Constructs real-valued frequency-domain data matrix 
%\Phi_L(V_{[0,M-1]}) = [F_L(V_{[0,M-1]}) F_L^*(V_{[1,M-1]})].
% Use this function for the single experiment/data set case.
% 
% INPUTS:
%   L : Depth
%   V : Frequency-domain data 
%   w : Frequency grid
% 
% OUTPUTS:
%   Phi : Frequency-domain data matrix \Phi_L(V_{[0,M-1]})
% 
% Date: October 28, 2025
% Author: T.J. Meijer - Eindhoven University of Technology
% Contact: t.j.meijer@tue.nl

% See the LICENSE file in the project root for full license information.

%% Construct \Phi_L(V_{[0,M-1]})
F = F_matrix(L,V,w);
Phi = [F conj(F(:,2:end))];

end