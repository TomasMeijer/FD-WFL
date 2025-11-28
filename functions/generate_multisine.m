function u = generate_multisine(w,P)
%GENERATE_MULTISINE Generates a random phase multisine signal
%
% INPUTS:
%   w : Excited frequencies
%   P : Number of periods

M = length(w); % Number of frequencies
phi = [0; 2*pi*rand(M-1,1)]; % Random phases (phi_0 = 0)
k = 0:2*M*P-1; % Bin vector

u = sum(cos(w.'*k+phi),1);

end