function Psi_real = Psi_matrix_real(L,V,w,E)
%PSI_MATRIX Constructs frequency-domain data matrix 
%\Psi_L(\{V^e_{[0,M-1]}\}_{e\in\mathcal{E}})T_{\Re} = [real(\mathcal{F}_L(\{V^e_{[0,M-1]}\}_{e\in\mathcal{E}})) imag(\mathcal{F}_L^*(\{V^e_{[0,M-1]}\}_{e\in\mathcal{E}}))].
% Use this function for the multi experiment/data set case.
% 
% INPUTS:
%   L : Depth
%   V : Frequency-domain data 
%   w : Frequency grid
%   E : Number of experiments/data sets
% 
% OUTPUTS:
%   Psi : Frequency-domain data matrix \Psi_L(\{V^e_{[0,M-1]}\}_{e\in\mathcal{E}})T_{\Re}
% 
% Date: October 28, 2025
% Author: T.J. Meijer - Eindhoven University of Technology
% Contact: t.j.meijer@tue.nl

% See the LICENSE file in the project root for full license information.

%% Construct \Psi_L(\{V^e_{[0,M-1]}\}_{e\in\mathcal{E}})T_{\Re}
Fcal = Fcal_matrix(L,V,w,E);
Psi_real = [real(Fcal) imag(Fcal(:,E+1:end))];

end