function Psi = Psi_matrix(L,V,w,E)
%PHI_MATRIX Constructs frequency-domain data matrix 
%\Psi_L(\{V^e_{[0,M-1]}\}_{e\in\mathcal{E}}) = [\mathcal{F}_L(V_{[0,M-1]}) \mathcal{F}_L^*(V_{[1,M-1]})].
% Use this function for the multi experiment/data set case.
% 
% INPUTS:
%   L : Depth
%   V : Frequency-domain data 
%   w : Frequency grid
%   E : Number of experiments/data sets
% 
% OUTPUTS:
%   Psi : Frequency-domain data matrix \Psi_L(\{V^e_{[0,M-1]}\}_{e\in\mathcal{E}})
% 
% Date: October 28, 2025
% Author: T.J. Meijer - Eindhoven University of Technology
% Contact: t.j.meijer@tue.nl

% See the LICENSE file in the project root for full license information.

%% Construct \Psi_L(\{V^e_{[0,M-1]}\}_{e\in\mathcal{E}})
Fcal = Fcal_matrix(L,V,w,E);
Psi = [Fcal conj(Fcal(:,E+1:end))];

end