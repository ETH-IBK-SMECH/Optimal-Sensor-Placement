function  [sensorDOF,varargout] = optimalSensorPlacement(un,wn,sensorNo,varargin)

% OPTIMALSENSORPLACEMENT perform optimal sensor placement
%
%   sensorDOF = optimalSensorPlacement(un,wn,sensorNo)
%
%   returns the degrees of freedom that must be instrumented, in 
%   order to measure the mode shapes un of a structure, by using 
%   sensorNo number of sensors. wn must be vector of the natural 
%   frequencies (in rad/s) that correspond to the mode shapes un. 
%
%   The above syntax implements the effective independence method
%   for the sensor setup. You can select another method using the 
%   syntax 
%
%   sensorDOF = optimalSensorPlacement(un,wn,sensorNo,'method',flg)
%
%   where flg specifies the method to be implemented. Select among:
%
%   'efi'     , for the effective independence (EFI, default)
%   'dpr'     , for the driving point residue (DPR)
%   'efi-dpr' , for the combined EFI-DPR
%   'mke'     , for the maximum kinetic energy (MKE)
%   'iei',    , for the information enthropy index(IEI)
%
%   sensorDOF = optimalSensorPlacement(un,wn,sensorNo,'method','iei','sigma',S)
%
%   implements a prediciton-error version of the OSP problem, where
%   S is a nondiagonal, positive-definite and symmetric correlation 
%   matrix.
%
%   [sensorDOF,iei] = optimalSensorPlacement(...
%
%   further returns a two column matrix, each row of which contains 
%   the information entropy index for a sensor placement at a given
%   number of DOFs. 
%
%   optimalSensorPlacement(...
%
%   plots the iei values against the number of sensors.
%
%   See also sap2Modal, sap2osp.
%

%
% Author: L. Claude & V. Dertimanis
% 1st Ed: 25-06-2014
% Last Update: 31-10-2016
% ETH Zurich 
% Institute of Structural Engineering 
% Chair of Structural Mechanics
%

% check inputs
if size(un,2) > sensorNo
    warning('smech:optimalSensorPlacement','Sensors are too few for measuring this number of mode shapes.')
end
if size(wn,1) ~= 1; wn = wn.';  end % we need wn to be row vector

% assign default options
defaultopt = struct('method','efi','sigma',eye(size(un,1)));
if nargin > 3
    % check how many varargins are present and assign data
    nvarargin = nargin-3;
    if any(nvarargin == [2 4]) == 0
        error('smech:optimalSensorPlacement','Please enter function parameters PROPERTY/VALUE pairs.')
    end
    for INDEX1 = 1:2:nvarargin        
        if ischar(varargin{INDEX1}) == 0
            error('smech:optimalSensorPlacement','Please enter function parameters PROPERTY/VALUE pairs.')
        end        
        switch varargin{INDEX1}             
            case 'method'
                if ischar(varargin{INDEX1+1}) == 0
                    error('smech:sap2osp','Invalid OSP method.')
                end
                defaultopt.method = varargin{INDEX1+1};
            case 'sigma'
                defaultopt.sigma = varargin{INDEX1+1};
        end        
    end
end

% process inputs according to method
PHI = un;
sensorDOF = (1:size(PHI,1))';
INDEX1 = 1;
L = eye(size(PHI,1));
S = defaultopt.sigma; 
R = chol(L*S*L'); 
iS = R\(R'\eye(size(L*S*L')));
FIM = (L*PHI)'*iS*(L*PHI);
refFIM = FIM;
informationEntropyIndex(:,1) = sqrt(det(refFIM))*ones(size(PHI,1)-sensorNo+1,1); % stores the IEI criterion
informationEntropyIndex(1,1) = 1;
informationEntropyIndex(:,2) = size(PHI,1):-1:sensorNo;                          % stores the sensorNo per iteration 
    
while size(PHI,1) > sensorNo    
    
    switch defaultopt.method
        
        case 'efi'
            [psi,lam] = eig(FIM,'nobalance');
            % effective independence distribution vector
            ilam = diag(1./diag(lam));
            QQ = (PHI*psi).*(PHI*psi)*ilam;
            informationVector = sum(QQ,2);
            
        case 'dpr'
            DPR = (PHI).*(PHI*diag(wn));
            minVal = min(DPR,[],2);
            informationVector = sum(DPR,2).*minVal;
            
        case 'efi-dpr'
            [psi,lam] = eig(FIM,'nobalance');
            % effective independence distribution vector
            ilam = diag(1./diag(lam));
            QQ = (PHI*psi).*(PHI*psi)*ilam;
            ED = sum(QQ,2);
            % DPR weights
            WN = repmat(wn,size(PHI,1),1);
            DPR = sum((PHI.^2)./WN,2);
            % scaled ED
            informationVector = ED.*DPR;
            
        case 'mke'
            MKE = PHI.*PHI;
            minVal = min(MKE,[],2);
            informationVector = sum(MKE,2).*minVal;
            
        case 'iei' 
            informationVector = zeros(1,size(sensorDOF,1));
            for INDEX2 = 1:1:size(sensorDOF,1)
                L = eye(size(PHI,1));
                % remove the sensor at position i
                L(INDEX2,:)=[];
                % calculate the iei-based FIM
                iR = chol(L*S*L'); 
                iiS = iR\(iR'\eye(size(L*S*L')));
                iFIM = (L*PHI)'*iiS*(L*PHI);
                informationVector(1,INDEX2) = 1/2*size(PHI,2)*log(2*pi())-1/2*log(det(iFIM));
            end
            
        otherwise
            error('smech:optimalSensorPlacement','Unknown OSP method.')
            
    end
    
    % get DOF with the minimum DPR-ED entry
    [~,I] = min(informationVector);
    % delete DOF, update data and return to top
    sensorDOF = [sensorDOF(1:I-1,1);sensorDOF(I+1:end,1)];
    PHI = [PHI(1:I-1,:);PHI(I+1:end,:)];
    L = eye(size(PHI,1));
    S = [S(1:I-1,:);S(I+1:end,:)];     
    S = [S(:,1:I-1),S(:,I+1:end)];    
    R = chol(L*S*L'); 
    iS = R\(R'\eye(size(L*S*L')));
    FIM = (L*PHI)'*iS*(L*PHI);
    informationEntropyIndex(INDEX1+1,1) = informationEntropyIndex(INDEX1+1,1)/sqrt(det(FIM));
    INDEX1 = INDEX1 + 1;
    
end

if nargout == 0
    % plot IEI for the implemented method
    bar(informationEntropyIndex(:,2),log10(informationEntropyIndex(:,1)),'k')
    hold on;
    bar(informationEntropyIndex(end,2),log10(informationEntropyIndex(end,1)),'r')
    xlabel('sensor number')
    ylabel('log_{10}(IEI)')
    set(gca,'xdir','reverse')
    axis tight
end

varargout(1) = {informationEntropyIndex};





