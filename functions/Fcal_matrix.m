function Fcal = Fcal_matrix(L,V,w,E)
%FCAL_MATRIX Constructs frequency-domain data matrix \mathcal{F}_L(\{V^e_{[0,M-1]}\}_{e\in\mathcal{E}}).
% Use this function for the single experiment/data set case.
% 
% INPUTS:
%   L : Depth
%   V : Frequency-domain data 
%   w : Frequency grid
%   E : Number of experiments/data sets
% 
% OUTPUTS:
%   Fcal : Frequency-domain data matrix \mathcal{F}_L(\{V^e_{[0,M-1]}\}_{e\in\mathcal{E}})
% 
% Date: October 28, 2025
% Author: T.J. Meijer - Eindhoven University of Technology
% Contact: t.j.meijer@tue.nl

% See the LICENSE file in the project root for full license information.

%% Prerequisites
W_L = @(w,L) exp(1i*w*[0:L-1].'); % Vector of increasing phase shifts
nv = size(V,1); % Dimension of V_k
M = length(w); % Number of frequencies in w (excluding negative frequencies)

%% Construct \mathcal{F}_L(\{V^e_{[0,M-1]}\}_{e\in\mathcal{E}})
Fcal = zeros(L*nv,M*E);
for m = 1:M
    for e = 1:E
        Fcal(:,(m-1)*E+e) = kron(W_L(w(m),L),V(:,m,e));
    end
end

end