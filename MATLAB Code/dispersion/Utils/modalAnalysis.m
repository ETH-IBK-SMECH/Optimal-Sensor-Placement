function [wn,zn,un] = modalAnalysis(M,C,K)

%
% modalAnalysis
%
% performs modal analysis of structural systems.
%
% [wn,zn,un] = modalAnalysis(M,C,K)
%
% Given the mass (M), stiffness (K) and damping (C) matrices, 
% modalAnalysis returns the natural frequncies wn (in rad/s), 
% the damping ratios zn and the mode shapes un.
%
% If the function is executed without output arguments, prints
% the results on workspace.
%
% See also dispersionAnalysis.
%
% Reference page in Help browser:
% <a href="matlab: web([docroot '/toolbox/mdac/funref/modalAnalysis.html'],'-helpbrowser')">doc modalAnalysis</a>
%

%
% Author: V. Ntertimanis
% 1st Ed: 28-10-2006
% Last Update: 31-05-2011
% National Technical University of Athens
% School of Mechanical Engineering
% Department of Mechanical Design & Automatic Control
% Copyright 1995-2011 V.K. Ntertimanis
%

% check damping
[P,L,C_,dampingStatus] = checkDamping(M,C,K);
if dampingStatus == 1     
    % undamped, or proportionally damped system    
    wn = sqrt(diag(L)); 
    zn = diag(C_)./(2*wn);
    un = M^(-1/2)*P;
    % store in ascending frequency    
    [wn,ind] = sort(wn);
    zn = zn(ind);
    un = un(:,ind');
else
    % (non proportionally) damped system    
    % solve a 2n eigenvalue problem
    A = [C M;M zeros(size(M))];
    B = [K zeros(size(M));zeros(size(M)) -M];
    [P,L] = eig(-B,A); 
    % obtain natural frequencies and damping ratios
    ReL = (L+L')/2;
    ImL = (L'-L)*1i/2;
    wn = (ReL^2+ImL^2).^.5;
    zn = -ReL./wn;
    % reform
    wn = diag(wn);wn = wn(1:2:end);
    zn = diag(zn);zn = zn(1:2:end);
    un = P(1:size(M,1),:);
end

if nargout == 0 
    
    form = '%7.3f';
    modeNoStruc = dprint((1:length(wn))','Mode No.','%1.0f');
    % reform modal data for display
    if dampingStatus == 1
        % print eigenvalues
        eigenValuesStruc = dprint(wn.^2,'Eigenvalue',form);
        % normalize eigenvectors
        for k = 1:size(un,2)
            [v,ind] = max(abs(un(:,k)));
            un(:,k) = un(:,k)/(sign(un(ind,k))*v);
        end
        % print eigenvectors
        s = dprint(un(:,1),'Mode 1',form);
        for k = 2:length(wn)
            modeNo = ['Mode ',num2str(k)];
            modePrint = dprint(un(:,k),modeNo,form);
            s = [s modePrint];
        end
    else
        % complex eigenvalues / eigenvectors
        eigVals = diag(L);
        eigReal = num2str(real(eigVals(1:2:end)),form);
        eigImag = num2str(abs(imag(eigVals(1:2:end))),form);
        eigStruc = blanks(size(eigReal,2)+size(eigImag,2)+6);       
        eigVecs = un(:,1:2:end);
        eigVecCellStruc = cell(1,size(eigReal,1));
        for k = 1:size(eigReal,1)
            % create structured vector of eigenvalues
            eigStruc(k,:) = [eigReal(k,:),' +/- ',eigImag(k,:),'i'];
            % create cell structured matrix of eigenvectors
            eigVecReal = num2str(real(eigVecs(:,k)),form);
            eigVecImag = num2str(abs(imag(eigVecs(:,k))),form);
            eigVecStruc = blanks(size(eigVecReal,2)+size(eigVecImag,2)+6);
            for kk = 1:size(eigReal,1)
                eigVecStruc(kk,:) = [eigVecReal(kk,:),' +/- ',eigVecImag(kk,:),'i'];
            end
            eigVecCellStruc(1,k) = {eigVecStruc};
        end
        % print eigenvalues
        eigenValuesStruc = dprint(eigStruc,'Eigenvalues',form);
        % print eigenvectors
        s = dprint(eigVecCellStruc{1,1},'Mode 1',form);
        for k = 2:length(wn)
            modeNo = ['Mode ',num2str(k)];
            modePrint = dprint(eigVecCellStruc{1,k},modeNo,form);
            s = [s modePrint];
        end        
    end
    wnStruc = dprint(wn/(2*pi),'Natural Frequency (Hz)',form);
    znStruc = dprint(100*zn,'Damping Ratio (%)',form);    
    disp(' ');
    disp(' ');
    disp('Natural frequencies and damping ratios')
    disp('---------------------------------------')
    disp([modeNoStruc wnStruc znStruc eigenValuesStruc]);
    disp(' ');
    disp('Mode Shapes')
    disp('---------------------------------------')
    disp(s)
    disp(' ')    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------
function [P,L,C_,dampingStatus] = checkDamping(M,C,K)

%
% calculates the symmetric eigenvalue
% problem and checks proportionality.
%


% calculate mass normalized stiffness
K_ = M^(-1/2)*K*M^(-1/2);
K_ = (K_ + K_.')/2; % to avoid numerical errors;
% solve the symmetric eigenvalue problem
[P,L] = eig(K_);
% attempt to diagonalize C
C_ = P'*M^(-1/2)*C*M^(-1/2)*P;
if max(max(abs(triu(C_,1)))) < sqrt(eps)
    dampingStatus = 1;
else
    dampingStatus = 0;
end
%--------------------------------------



