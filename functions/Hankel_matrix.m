function H = Hankel_matrix(L,v)
%HANKEL_MATRIX Constructs a Hankel matrix of depth L from the sequence v_{[0,N-1]}
%H_L(v_{[0,N-1]}) = [v_{[0,L-1]} v_{[1,L]} ... v_{[N-L+1,N-1]}].
% 
% INPUTS:
%   L : Depth
%   V : Frequency-domain data 
%   E : Number of experiments/data sets
% 
% OUTPUTS:
%   H : Hankel matrix H_L(v_{[0,N-1]})
% 
% Date: November 27, 2025
% Author: T.J. Meijer - Eindhoven University of Technology
% Contact: t.j.meijer@tue.nl

% See the LICENSE file in the project root for full license information.

[nv,N,E] = size(v);

%% Construct H_L(v_{[0,N-1]})
H = zeros(L*nv,E*(N-L+1));
for e = 1:E
    for i = 1:N-L+1
        v_seg = v(:,(0:L-1)+i,e);
        H(:,(e-1)*(N-L+1)+i) = v_seg(:);
    end
end
end

