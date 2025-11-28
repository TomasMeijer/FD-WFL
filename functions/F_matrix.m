function F = F_matrix(L,V,w)
%F_MATRIX Constructs frequency-domain data matrix F_L(V_{[0,M-1]}).
% Use this function for the single experiment/data set case.
% 
% INPUTS:
%   L : Depth
%   V : Frequency-domain data 
%   w : Frequency grid
% 
% OUTPUTS:
%   F : Frequency-domain data matrix F_L(V_{[0,M-1]})
% 
% Date: October 28, 2025
% Author: T.J. Meijer - Eindhoven University of Technology
% Contact: t.j.meijer@tue.nl

% See the LICENSE file in the project root for full license information.

%% Prerequisites
W_L = @(w,L) exp(1i*w*[0:L-1].'); % Vector of increasing phase shifts
nv = size(V,1); % Dimension of V_k
M = length(w); % Number of frequencies in w (excluding negative frequencies)

%% Construct F_L(V_{[0,M-1]})
F = zeros(L*nv,M);
for m = 1:M
    F(:,m) = kron(W_L(w(m),L),V(:,m));
end

end