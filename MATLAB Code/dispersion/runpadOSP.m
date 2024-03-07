% This script carries out a disperions analysis for a FEM model and
% computes optimal sensor configurations using several methods (MKE, EFI,
% IEI)


% Authors: V. Dertimanis & C. Leyder  
% 1st Ed: 25-06-2014
% Last Update: 04-12-2018
% ETH Zurich 
% Institute of Structural Engineering 
% Chair of Structural Mechanics

addpath('./Utils')

% %% SECTION 1: READ SAP2000 MATRICES (OR FROM ANY OTHER TYPE OF FEM MODEL)
% %
% % This section:
% %
% % 1. reads the excel file where the SAP2000 tx* files are stored in
% % sheets (=element, mass and stiffness matrix from the FEM model)
% % 2. creates the mass and stiffness matrix of the SAP2000 2D frame
% % 3. performs static condensation and removes the massless DOFs
% % 4. creates a proportional damping matrix for the reduced system
% % 4. creates the input matrix
% % 5. saves the reduced matrices
% %
% 
% clear
% close all
% clc
% 
% % excitation: select among the following
% e1 = 22;  % corresponds to equation No. for the U3 direction of node ~25 (left beam)
% e2 = 2;   % corresponds to equation No. for the U3 direction of node ~38 (center beam)
% e3 = 156; % corresponds to equation No. for the U3 direction of node ~50 (right beam)
% e4 = 57;  % corresponds to equation No. for the U1 direction of node ~3  (column 1)
% e5 = 54;  % corresponds to equation No. for the U1 direction of node ~13 (column 2)
% e6 = 147; % corresponds to equation No. for the U1 direction of node ~18 (column 3)
% e7 = 150; % corresponds to equation No. for the U1 direction of node ~8  (column 4)
% e = [e2 e4]; % adjust accordingly
% 
% % provide damping
% zn = 0.015;
% 
% % mass data
% txm = xlsread('sap2000TXFiles.xlsx','txm','A1:C184');
% % stiffness data
% txk = xlsread('sap2000TXFiles.xlsx','txk','A1:C907');
% [~,~,txe] = xlsread('sap2000TXFiles.xlsx','txe','A1:F64');
% U1_eq=zeros(size(txe,1),1); U3_eq=U1_eq; R2_eq=U1_eq;
% 
% 
% %read equation numbers and DOF data
% DOF_labels=cell(length(U1_eq),1);
% 
% for i=1:length(U1_eq)
% U1_eq(i,1)=txe{i,2};
% U3_eq(i,1)=txe{i,4};
% R2_eq(i,1)=txe{i,6};
% DOF_labels{i,1}=txe{i,1}; %these DOF labels are already different from the DOF labels used by CL
% end
% % DOF_numbering_CL(:,1)=[1,2,7,8,3,4,5,6,9:64];
% 
% DOF_dir=zeros(length(U1_eq),1);
% 
% DOF_K = cell(192-8,1);   %lenght of CL = 192 - 4 supports * 2 !!! %Create a list with which DOF is where in the stiffness matrix
% for IND = 1:length(U1_eq)
%     if U1_eq(IND,1)~=0 
%     DOF_K(U1_eq(IND,1)) = DOF_labels(IND,1);
%     DOF_dir(U1_eq(IND,1),1)=U1_eq(IND,1);
%     end
%     if U3_eq(IND,1)~=0 
%     DOF_K(U3_eq(IND,1)) = DOF_labels(IND,1);
%     DOF_dir(U3_eq(IND,1),2)=U3_eq(IND,1);
%     end
%     if R2_eq(IND,1)~=0
%     DOF_K(R2_eq(IND,1)) = DOF_labels(IND,1);
%     DOF_dir(R2_eq(IND,1),3)=R2_eq(IND,1);
%     end
% end
% 
% 
% % mass and stiffness matrices
% M = zeros(length(txm)); K = M;
% for IND = 1:length(txm)
%     M(txm(IND,1),txm(IND,2)) = txm(IND,3);
% end
% for IND = 1:length(txk)
%     K(txk(IND,1),txk(IND,2)) = txk(IND,3);
% end
% K = K + tril(K,-1).'; %completes the triangle
% 
% % 
% 
% % permute M and K in order to bring first the massless DOFs
% s = find(diag(M) == 0);
% p = find(diag(M) ~= 0);
% permutationVector = [s;p];
% M = M(permutationVector,permutationVector);
% K = K(permutationVector,permutationVector);
% DOF_K=DOF_K(permutationVector);
% DOF_dir=DOF_dir(p,:);
% 
% DOF_K =DOF_K(length(s)+1:end);
% 
% 
% % static condensation
% Kss = K(1:length(s),1:length(s));
% Ksp = K(1:length(s),length(s)+1:end);
% T_ = -Kss\Ksp;
% T = [T_;eye(length(p))];
% K = T'*K*T;
% % K = (K + K.')/2; % to avoid numerical errors;
% M = M(length(s)+1:end,length(s)+1:end);
% C = getDamping(M,K,zn*ones(length(p),1));
% % C = (C + C.')/2; % to avoid numerical errors;
% % input matrix
% P = zeros(length(p),length(e));
% for IND = 1:length(e)
%     fi = find(p(:) == e(IND));
%     P(fi,IND) = 1;
% end
% 
% save sapMatrices M C K P p DOF_K DOF_dir



% %% SECTION 2: CALCULATE ENERGY FOR THE TOTAL DOFs (ACCELERATION SENSORS)
% %
% % This section
% %
% % 1. loads the structural matrices of the 2D frame
% % 2. creates the state-space model of the structure
% % 3. performs dispersion analysis for the total DOF output
% % 4. saves the trace of the zero-lag covariance matrix of the output
% %
% 
% clear
% close all
% clc
% 
% senstype='disp';  %distinguish between acceleration, velocity and displacement sensors (different disperion metrics)
% 
% % load frame
% load sapMatrices
% n = size(M,1);
% m = size(P,2);
% 
% % create discrete-time state-space
% Ts = 1/(5*max(sqrt(eig(K,M))/(2*pi)));
% Ac = [zeros(n) eye(n);-M\K -M\C];
% Bc = [zeros(n,m);M\P];
% %displacement
% if strcmp(senstype,'disp')
% Cc = [eye(n) zeros(n)]; 
% Dc = 0;                 
% end
% % %velocity
% if strcmp(senstype,'vel')
% Cc = [zeros(n) eye(n)];  
% Dc = 0;                  
% end
% % %acceleration
% if strcmp(senstype,'accel')
% Cc=[-M\K -M\C];
% Dc=M\P;
% end
% 
% if strcmp(senstype,'accel')
% totalStateSpace_acc = ss(Ac,Bc,Cc,Dc);
% end
% if strcmp(senstype,'vel')
% totalStateSpace_vel = ss(Ac,Bc,Cc,Dc);
% end
% if strcmp(senstype,'disp')
% totalStateSpace_disp = ss(Ac,Bc,Cc,Dc);
% end
% 
% % totalStateSpace = c2d(totalStateSpace,Ts); %discrete time
% 
% % dispersion analysis
% sigmaFF = eye(m);
% if strcmp(senstype,'accel')
% [~,~,Dacc,~] = dispersionAnalysis(totalStateSpace_acc,sigmaFF);  %run without output to get table
% 
% % energy metric all dofs measured
% totalEnergyMetric_acc = diag(real(sum(Dacc,3)));
% 
% % save metric
% save totalEnergyMetric_acc totalEnergyMetric_acc
% save totalStateSpace_acc totalStateSpace_acc sigmaFF
% end
% if strcmp(senstype,'vel')
% [~,~,Dvel,~] = dispersionAnalysis(totalStateSpace_vel,sigmaFF);  %run without output to get table
% 
% % energy metric all dofs measured
% totalEnergyMetric_vel = diag(real(sum(Dvel,3)));
% 
% % save metric
% save totalEnergyMetric_vel totalEnergyMetric_vel
% save totalStateSpace_vel totalStateSpace_vel sigmaFF
% end
% if strcmp(senstype,'disp')
% [~,~,Ddisp,~] = dispersionAnalysis(totalStateSpace_disp,sigmaFF);  %run without output to get table
% 
% % energy metric all dofs measured
% totalEnergyMetric_disp = diag(real(sum(Ddisp,3)));
% 
% % save metric
% save totalEnergyMetric_disp totalEnergyMetric_disp
% save totalStateSpace_disp totalStateSpace_disp sigmaFF
% end



%% SECTION 3: OSP (so far, for acceleration sensors only)
% 
% This section:
%
% 1. solves the optimal sensor placement problem for various sensor numbers
% 2. for each sensor number, calculates the dispersion-based energy metric
% 3. plots the results for the each OSP method employed
%

clear
close all
clc

senstype='disp';  %'disp','vel' specify different sensor type



% set sensor range
sensorRange = 6:28;

% set modes for osp
ospModes = 1:5; 

% load structural matrices
load sapMatrices
n = size(M,1); % dof

if strcmp(senstype,'acc')
% load total state-space model
load totalStateSpace_acc

% load dispersion-based metric (for acceleration sensors)
load totalEnergyMetric_acc

elseif strcmp(senstype,'vel')

load totalStateSpace_vel
load totalEnergyMetric_vel   
    
elseif strcmp(senstype,'disp')
    
load totalStateSpace_disp
load totalEnergyMetric_disp       
    
end

%load covariance matrix for IEI method and convert to current DOF numbering
cd ../
load Sigma
cd 'dispersion'

%plotting specifications
grey=[0.65 0.65 0.65];

% %Build up new Sigma matrix (new numbering of DOFs)
perm_j=zeros(112,1);

%remove tabs from DOF_K
for i=[1:44,47:56,59:78,81:112,115:120]
DOF_K{i,1} = regexprep(DOF_K{i,1}, '\s', ''); %#ok<SAGROW>
end

for i=[1:44,47:56,59:78,81:112,115:120]
perm_j(i,1)=find(strcmp(DOF_K(i),Coord.Name)); 
end
perm_j(45:46,1)=[2,2];
perm_j(57:58,1)=[4,4];
perm_j(79:80,1)=[8,8];
perm_j(113:114,1)=[6,6];

perm_dof=zeros(120,1);
for i=1:1:length(perm_j)
%x-direction
if DOF_dir(i,1)~=0
perm_dof(i)=(perm_j(i)-1)*3+1;
end
%y-direction
if DOF_dir(i,2)~=0
perm_dof(i)=(perm_j(i)-1)*3+2;
end
end

sigmaN=sigma(perm_dof,perm_dof);


% calculate modes
[wn,~,un] = modalAnalysis(M,0,K);

for i=1:1:length(wn)  %norm all mode shapes to 1 (every mode shape is considered with equal importance) 
    [mval, ind]=max(abs(un(:,i)));
    un(:,i)=un(:,i)*sign(un(ind,i))/mval;
end

fn=wn./(2*pi);

% loop over sensorNo
mkeMetric = zeros(1,length(sensorRange));
efiMetric = zeros(1,length(sensorRange));
ieiMetric = zeros(1,length(sensorRange));

R=struct;
R_vd=struct;

for IND = 1:length(sensorRange)
    
    % display current sensor
    disp(['Sensor No.: ',num2str(sensorRange(IND))]);
    
    % solve osp
    R(sensorRange(IND)).ospMKE = optimalSensorPlacement(un(:,ospModes),wn(ospModes),sensorRange(IND),'method','mke');
    R(sensorRange(IND)).ospEFI = optimalSensorPlacement(un(:,ospModes),wn(ospModes),sensorRange(IND),'method','efi','sigma',sigmaN);
    R(sensorRange(IND)).ospIEI = optimalSensorPlacement(un(:,ospModes),wn(ospModes),sensorRange(IND),'method','iei','sigma',sigmaN);
    %note: use perm_dof(ospMKE) to convert to DOF number from main code
%     [sensorDOF,FIM,PHI,H,infoVec,pos,sigma,Accuracy,un,EFI_pareto] = OSP_COR_PE('iei',un(:,ospModes),wn(ospModes),sensorRange(IND),sigmaN);
% both functions deliver the same results, the difference lies in the
% input!!!

    %permute to original coordinate system
    R_vd(sensorRange(IND)).ospMKE=perm_dof(R(sensorRange(IND)).ospMKE);
    R_vd(sensorRange(IND)).ospIEI=perm_dof(R(sensorRange(IND)).ospIEI);
    R_vd(sensorRange(IND)).ospEFI=perm_dof(R(sensorRange(IND)).ospEFI); 

if strcmp(senstype,'accel')
    mkeMetric(sensorRange(IND)) = 100*sum(totalEnergyMetric_acc(R(sensorRange(IND)).ospMKE))/sum(totalEnergyMetric_acc);
    efiMetric(sensorRange(IND)) = 100*sum(totalEnergyMetric_acc(R(sensorRange(IND)).ospEFI))/sum(totalEnergyMetric_acc);
    ieiMetric(sensorRange(IND)) = 100*sum(totalEnergyMetric_acc(R(sensorRange(IND)).ospIEI))/sum(totalEnergyMetric_acc);
elseif strcmp(senstype,'vel')
    mkeMetric(sensorRange(IND)) = 100*sum(totalEnergyMetric_vel(R(sensorRange(IND)).ospMKE))/sum(totalEnergyMetric_vel);
    efiMetric(sensorRange(IND)) = 100*sum(totalEnergyMetric_vel(R(sensorRange(IND)).ospEFI))/sum(totalEnergyMetric_vel);
    ieiMetric(sensorRange(IND)) = 100*sum(totalEnergyMetric_vel(R(sensorRange(IND)).ospIEI))/sum(totalEnergyMetric_vel);
elseif strcmp(senstype,'disp')
    mkeMetric(sensorRange(IND)) = 100*sum(totalEnergyMetric_disp(R(sensorRange(IND)).ospMKE))/sum(totalEnergyMetric_disp);
    efiMetric(sensorRange(IND)) = 100*sum(totalEnergyMetric_disp(R(sensorRange(IND)).ospEFI))/sum(totalEnergyMetric_disp);
    ieiMetric(sensorRange(IND)) = 100*sum(totalEnergyMetric_disp(R(sensorRange(IND)).ospIEI))/sum(totalEnergyMetric_disp);
end

        
end

% plot the mke, efi and iei energy metrics (separate plots)
figure
subplot(311)
bar(sensorRange,nonzeros(mkeMetric)), axis tight
ylabel('MKE energy ratio (%)')
subplot(312)
bar(sensorRange,nonzeros(efiMetric)), axis tight
ylabel('EFI energy ratio (%)')
subplot(313)
bar(sensorRange,nonzeros(ieiMetric)), axis tight
xlabel('sensor number')
ylabel('Energy ratio (%)')


% plot the mke, efi and iei energy metrics (same plot)
figure
plot(sensorRange,nonzeros(mkeMetric),'k*-') 
axis tight;
hold on; grid on;
plot(sensorRange,nonzeros(efiMetric),'*-','Color',grey)
plot(sensorRange,nonzeros(ieiMetric),'ko-')
xlabel('sensor number')
ylabel('Energy ratio (%)')
legend('MKE','EFI','IEI')

cd ../
save('mac calculation\OSP_vd.mat','R_vd','sensorRange')
cd 'dispersion'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For random Configurations

cd ..\
load 'random configurations/RandomConfigs.mat' randind
%randind: random samples with unique values (sampling without replacement!!!) - generated by Generate_randomConfig.m (DOF numbering according to Comparison_OSP.m)
%CONVERT TO DOF NUMBERING FROM THIS FUNCTION!!!
%randind(DOFs,Samples,nSensors)

[~,nsamples,~]=size(randind); 
randMetric = zeros(max(sensorRange),nsamples);

for nsensors=min(sensorRange):max(sensorRange)
  for i=1:nsamples
      clear prior post
      prior=randind(1:nsensors,i,nsensors-min(sensorRange)+1); %DOF numbering prior to conversion
      post=zeros(length(prior),1);
      for j=1:1:length(prior)
          if ~isempty(find(perm_dof==prior(j),1))
          post(j)=find(perm_dof==prior(j)); %DOF numbering posterior to conversion (DOFs that are not found are support DOFs, their energy contribution is 0 anyway)
          end
      end
        post(post==0)=[];
        if strcmp(senstype,'accel')
        randMetric(nsensors,i) = 100*sum(totalEnergyMetric_acc(post))/sum(totalEnergyMetric_acc); 
        elseif strcmp(senstype,'vel')
        randMetric(nsensors,i) = 100*sum(totalEnergyMetric_vel(post))/sum(totalEnergyMetric_vel); 
        elseif strcmp(senstype,'disp')
        randMetric(nsensors,i) = 100*sum(totalEnergyMetric_disp(post))/sum(totalEnergyMetric_disp); 
        end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% For an Implemented Configuration (EXAMPLE)
% Jointno=[2,11,34,4,21,47,6,26,59,8,16];
% DofNox=(Jointno-1)*3+1;
% DofNoy=(Jointno-1)*3+2;
% Dofno=[DofNox,DofNoy];
% for i=1:1:length(Dofno)
%     DofnoVD(i)=find(perm_dof==Dofno(i));
% end
% 
% configmetric = 100*sum(totalEnergyMetric_acc(DofnoVD))/sum(totalEnergyMetric_acc);


%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% save Dispersion metrics 

if strcmp(senstype,'acc')
Metrics.DispAcc.rand=randMetric;

Metrics.DispAcc.MKE=mkeMetric;
Metrics.DispAcc.EFI=efiMetric;
Metrics.DispAcc.IEI=ieiMetric;

save('paretofront\MetricsDispAcc.mat','Metrics')

elseif strcmp(senstype,'vel')

Metrics.DispVel.rand=randMetric;

Metrics.DispVel.MKE=mkeMetric;
Metrics.DispVel.EFI=efiMetric;
Metrics.DispVel.IEI=ieiMetric;

save('paretofront\MetricsDispVel.mat','Metrics')
    
elseif strcmp(senstype,'disp')
    
Metrics.DispDisp.rand=randMetric;

Metrics.DispDisp.MKE=mkeMetric;
Metrics.DispDisp.EFI=efiMetric;
Metrics.DispDisp.IEI=ieiMetric;

save('paretofront\MetricsDispDisp.mat','Metrics')     
    
end


cd 'dispersion'

