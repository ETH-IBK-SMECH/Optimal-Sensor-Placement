function  [sensorDOF,FIM,PHI,H,infoVec,pos,sigma,Accuracy] = OSP_COR_PE_MultiAxis(flg,un,wn,sensorNo,sigma)

% optimalSensorPlacement for tri-axial sensors (Core-function)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function INPUTs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% flg specifies the method to be implemented. Select among:
%
% 'efi'     , for the effective independence (EFI) method
% 'dpr'     , for the driving point residue (DPR) method
% 'efi-dpr' , for the combined EFI-DPR method
% 'mke'     , for the maximum kinetic energy (MKE) method
% 'iei',    , for the information enthropy index method (IEI)
%
%Remark: all algorithms use the 'BSSP' algorithm,where first a sensor at
%every DOF is assumed, and then the sensors are sequentially removed
%
% un:       mode shape matrix
% 
% wn:       must be a row vector of the natural frequencies (in rad/s) 
%           that correspond to the mode shapes un.
%
% sensorNo: number of DOFs to be instrumented / number of Sensors
%
% sigma:    correlation matrix, calculated in Evaluation_OSP.m 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Function OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%sensorDOF:     index number of the DOFs to be monitored
%FIM :          FIsher Information matrix (reduced to the monitored DOFs)
%PHI:           mode shape matrix (reduced to the monitored DOFs)
%H:             information entropy matrix (reduced to the monitored DOFs)
%infoVec:       matrix with the information entropy at every iteration step
%pos:           matrix with the monitored DOFs at every iteration step
%sigma:         correlation matrix (reduced to the monitored DOFs)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: V. Ntertimanis (modified by C. Leyder)
% 1st Ed: 25-06-2014
% Last Update: 21-11-2018 (CL)
% ETH Zurich
% Copyright 1995-2018 V.K. Ntertimanis, C. Leyder

if size(un,2) > sensorNo
    warning('mdac:optimalSensorPlacement','Sensors are too few for measuring this number of mode shapes.')
end 
PHI = un;
sensorDOF = (1:size(PHI,1))';
L_in=eye(size(PHI,1));
FIMref=(L_in*PHI)'*((L_in*sigma*L_in')\(L_in*PHI));

% positivedefinite = all(eig(FIMref) > 0) %use only for checking

k=0;
while size(PHI,1) > sensorNo    
    k=k+1;
    % eigenvalue problem of the Fisher Information Matrix
    FIM = PHI'*(sigma\PHI);
%     positivedefinite = all(eig(FIM) > 0) %use only for checking
    switch flg
        
        case 'efi'
            [psi,lam] = eig(FIM,'nobalance');   %[V,D] = eig(A) returns diagonal matrix D of eigenvalues and matrix V whose columns are the corresponding right eigenvectors, so that A*V = V*D.
            %lam: diagonal matrix of eigenvalues
            %psi: eigenvectors
            % effective independence distribution vector
%             Orthocheck=psi'*psi; %used only for checking
            ilam = diag(1./diag(lam)); %inverse of the diagonal matrix of eigenvalues
            QQ = (PHI*psi).*(PHI*psi)*ilam;
            informationVector = sum(QQ,2);
            EFI=informationVector;
            
        case 'dpr'
            DPR = (PHI).*(PHI*diag(wn));
%             minVal = min(DPR,[],2);
            informationVector = sum(DPR,2);%.*minVal;

        case 'efi-dpr'
            [psi,lam] = eig(FIM,'balance');
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
%             minVal = min(MKE,[],2);
            informationVector = sum(MKE,2);%.*minVal;
           % The iteration is superfluous for the mke method,
            %the values remain the same. Sensor locations could be chosen
            %immediately from step 1! For structure with uniform mass
            %distribution mke and efi should deliver the same results (i.e.
            %the mke is then the fast method, since no iteration is
            %required). For the test case frame, the mass is not evenly distributed
            %per DOF (columns and beams).
            %Iteration is nonetheless included for consistency with the
            %other methods
            
            %add EFI for pareto front plotting
            [psi,lam] = eig(FIM,'nobalance');   
            ilam = diag(1./diag(lam)); 
            QQ = (PHI*psi).*(PHI*psi)*ilam;
            EFI = sum(QQ,2);
            
            
        case 'iei' %via the information entropy
            clear IEI L L_new informationVector FIM_BSSP i H
            L_new=eye(size(PHI,1));   
            j=0;
            
            
            %pre-allocate:
            H=zeros(size(sensorDOF,1),1);
            IEI=zeros(size(sensorDOF,1),1);
            informationVector=zeros(size(sensorDOF,1),1);
            
            for i=1:3:size(sensorDOF,1) %step-size=3, for tri-axial sensors!
                j=j+1;
                L=L_new;
                L(i:i+2,:)=[];   %remove the sensor at position i (corresponds to 3 DOFs) and calculate the FIM & IEI of that configuration
                FIM_BSSP=(L*PHI)'*((L*sigma*L')\(L*PHI));
                %FIM_BSSP=(L*PHI)'*(L*sigma*L')^(-1)*(L*PHI); %very slow
                H(i,1)=1/2*size(PHI,2)*log(2*pi())-1/2*log(det(FIM_BSSP)); %find the best configuration, where Sensor i is removed (peaks in H are sensors that should not be removed, because then the information entropy would increase), downpeaks are sensors that should be removed, because these configurations are good compared to the others)
                IEI(i,1)=sqrt(det(FIMref)/det(FIM_BSSP));  %delete configuration with max determinant 
                informationVector(i,1)=H(i,1);
            end
            
            %add EFI for pareto front plotting
            [psi,lam] = eig(FIM,'nobalance');   
            ilam = diag(1./diag(lam)); 
            QQ = (PHI*psi).*(PHI*psi)*ilam;
            EFI = sum(QQ,2);
    %
        otherwise
            error('mdac:optimalSensorPlacement','Incorrect specification for the FLG argument.')
            
    end
    % get DOF with the minimum informationVector entry
    [~,I] = min(nonzeros(informationVector));
    I=(I-1)*3+1;  %remove three rows at once, corresponding to the same sensor (tri-axial)
    % delete DOF and return to top
    sensorDOF = [sensorDOF(1:I-1,1);sensorDOF(I+3:end,1)];
    PHI = [PHI(1:I-1,:);PHI(I+3:end,:)];
    sigma = [sigma(1:I-1,:);sigma(I+3:end,:)];
    sigma = [sigma(:,1:I-1),sigma(:,I+3:end)];
    informationVector=[informationVector(1:I-1,1);informationVector(I+3:end,1)];
    EFI=[EFI(1:I-1,1);EFI(I+3:end,1)];
    infoVec(1:length(informationVector),k)=informationVector; %#ok<AGROW>
    pos(1:length(sensorDOF),k)=sensorDOF; %#ok<AGROW>
    Accuracy(1,k)=det(FIM);     %#ok<AGROW>
end
FIM = PHI'*(sigma\PHI);
% positivedefinite = all(eig(FIM) > 0)
H=1/2*size(PHI,2)*log(2*pi())-1/2*log(det(FIM));
IEI=sqrt(det(FIMref)/det(FIM));     %#ok<NASGU>



