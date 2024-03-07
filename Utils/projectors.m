function [G,varargout] = projectors(A)

% 
% PROJECTORS
%
% spectral decomposition.
%
% G = projectors(A)
%
% calculates the spectral projectors of a
% square matrix A, which has complete set 
% of eigenvalues (called the spectrum). G
% returns a [n x n x k] array when A is a
% [n x n] matrix.
%
% Any [n x n] matrix A that has a complete
% set of eigenvalues di, i=1,...,k, may be
% decomposed as (k<=n):
%
% A = d1*G1+d2*G2+...+dk*Gk
%
% [G,d]=projectors(A)
%
% additionally returns the k eigenvalues D 
% of the matrix A.
%
% See also EIG.
%
% Reference page in Help browser:
% <a href="matlab: web([docroot '/toolbox/mdac/funref/projectors.html'],'-helpbrowser')">doc projectors</a>
%

%
% Author: V. Ntertimanis
% 1st Ed: 29-12-2007
% Last Update: 30-04-2013
% National Technical University of Athens
% Copyright 1995-2013 V.K. Ntertimanis
%

[n,m] = size(A);
if n ~= m
    error('mdac:projectors','Only square matrices are allowed.')
end
[V,D] = eig(A,'nobalance');
iV = V\eye(n);
G = zeros(n,n,n);
for k = 1:n
    G(1:n,1:n,k) = V(:,k)*iV(k,:);        
end
% outputs
varargout(1) = {diag(D)};

% 
% [n,m] = size(A);
% if n ~= m
%     error('mdac:projectors','Only square matrices are allowed.')
% end
% d(1:n,1) = eig(A,'nobalance');
% ds = singleeig(d,n);
% G = zeros(n,n,length(ds));
% for i = 1:length(ds)
%     Prod1 = eye(n);
%     Prod2 = 1;
%     for j = 1:length(ds)
%         if j ~= i
%             Prod1 = Prod1*(A-ds(j)*eye(n));
%             Prod2 = Prod2*(ds(i)-ds(j));
%         end
%     end
%     G(1:n,1:n,i) = Prod1/Prod2;    
% end
% % outputs
% varargout(1) = {ds};
% 
% 
% 
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function ds = singleeig(d,n)
% 
% %
% % check whether multiple eigenvalues
% % exist and return a vector of single
% % ones.
% %
% 
% d = roundnumber(d,1e06);%n=3;
% ind = 1;
% k = 1;
% ds = zeros(1,1);
% while ind < n+1
%     ds(k,1) = d(ind);
%     a = find(d == d(ind));
%     ind = ind+length(a);
%     k = k+1;
% end
% 
% 
% 
