% Evaluates the IEI Metric for the random configurations, generated with 
% random configurations/Generate_randomConfig.m

% Part0. for random configurations
% Part1. for a specific configuration (e.g. the implemented configuration)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: C. Leyder
% Last Update: 5.12.2018
% ETH Zurich
% Copyright 2018 C. Leyder

clearvars
close all
clc


cd ../
currentpath=pwd;
addpath(strcat(currentpath,'\mac calculation'))
addpath(strcat(currentpath,'\random configurations'))
cd 'paretofront'


load 'RandomConfigs.mat' randind          %DOF numbers from the random configurations
load 'IEI_Metric.mat' IEI_EFI IEI_MKE IEI_IEI                   %IEI metrics from the configurations derived with the OSP methods (mke,efi,iei)
load 'ModalData_FE.mat' atyp Coord un wn                              %modal data from the FE model (created within the Evaluation_OSP.m function)
load(strcat(currentpath,'\Sigma.mat'),'sigma')                   %Prediction error correlation matrix, to calculate the FIM of the reference configuration


%% 0. Evaluate the IEI metric for random configurations
[~,nsamples,~]=size(randind); %number of random configurations
[nsmin,~]=size(nonzeros(randind(:,1,1))); %minimum number of sensors
[nsmax,~]=size(nonzeros(randind(:,1,end))); %maximum number of sensors
[~,nm]=size(un); %number of target modes



clear IEI_rand
IEI_rand=zeros(nsmax,nsamples);

Qref=un'*(sigma\un);
for nsensors=nsmin:nsmax
    for i=1:1:nsamples
        r_ind=randind(1:nsensors,i,nsensors-nsmin+1);

        sens_ind.x=[];
        sens_ind.z=[];
        % ignore rotational DOFs (no tilt sensors)
        for j=1:1:nsensors
            if rem(r_ind(j),3)==2
                sens_ind.z(i,j)=r_ind(j);
            elseif rem(r_ind(j),3)==1
                sens_ind.x(i,j)=r_ind(j);
            end 
        end
        sens_ind.x(sens_ind.x==0)=[];
        sens_ind.z(sens_ind.z==0)=[];   
        sens_ind.x=sens_ind.x';
        sens_ind.z=sens_ind.z';
        
        r_ind_xz=[sens_ind.x;sens_ind.z];
        
        PHI=un(r_ind_xz,:);
        sigred=sigma(r_ind_xz,r_ind_xz);
        ConfigQ=PHI'*(sigred\PHI);
        IEI_rand(nsensors,i)=abs(sqrt(det(Qref)/det(ConfigQ)));
        clear ConfigQ PHI r_ind sigred
    end

end

save('IEI_Metric.mat','IEI_rand','IEI_MKE','IEI_EFI','IEI_IEI')


% %% 1. Evaluate for implemented configurations
% % FULL Configuration
% sens_ind.x=[2,11,34,4,21,47,6,26,59,8,16]';
% sens_ind.z=[2,11,34,4,21,47,6,26,59,8,16]';
% 
%         
% r_ind_xz=[(sens_ind.x-1)*3+1;(sens_ind.z-1)*3+2];
%         
% PHI=un(r_ind_xz,:);
% sigred=sigma(r_ind_xz,r_ind_xz);
% ConfigQ=PHI'*(sigred\PHI);
% IEI_config_full=abs(sqrt(det(Qref)/det(ConfigQ)));
% clear ConfigQ PHI r_ind sigred
% 
% 
% % REDUCED Configuration
% sens_ind.x=[2,11,34,4,21,47,6,26,59,8,16]';
% sens_ind.z=[34,47,59]';
%         
% r_ind_xz=[(sens_ind.x-1)*3+1;(sens_ind.z-1)*3+2];
%         
% PHI=un(r_ind_xz,:);
% sigred=sigma(r_ind_xz,r_ind_xz);
% ConfigQ=PHI'*(sigred\PHI);
% IEI_config_red=abs(sqrt(det(Qref)/det(ConfigQ)));
% clear ConfigQ PHI r_ind sigred



