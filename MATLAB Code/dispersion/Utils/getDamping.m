function C = getDamping(M,K,zn)

%
% getDamping
%
% estimate an extended Rayleigh viscous damping matrix.
%
% C = getDamping(M,K,zn)
%
% returns a proportional viscous matrix C, which corresponds to 
% modal damping values specified by the vector zn. Calculation
% is taking place by using the extended Rayleigh damping method.
%
% See also modalAnalysis.
%
% Reference page in Help browser:
% <a href="matlab: web([docroot '/toolbox/mdac/funref/getDamping.html'],'-helpbrowser')">doc getDamping</a>
%

%
% Author: V. Ntertimanis
% 1st Ed: 15-03-2014
% Last Update: 02-04-2014
% ETH Zurich
% Copyright 1995-2014 V.K. Ntertimanis
%

% 0. checks and creation of power vector
if length(zn) ~= length(M)
    error('mdac:getDamping','Inconsistency of input arguments.')
end
% 1. extract modals
[wn,~,un] = modalAnalysis(M,0,K);
% 2. generalized mass matrix
M_ = un'*M*un;
% 3. D matrix
D = diag(2*zn.*wn./diag(M_));
C = M*un*D*un'*M;
