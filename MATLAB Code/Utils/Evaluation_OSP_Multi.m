function [Results,Qref]=Evaluation_OSP_Multi(flg,filename,atyp,tm,ns,lambda,el_sel,ploton)
%Evaluates the optimal sensor placement algorithms based on the input
%defined in Comparison_OSP.m, defines addtional input parameters, such as
%the correlation matrix sigma (for 3-AXIAL sensors)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: C. Leyder
% Last Update: 21.11.2018
% ETH Zurich
% Copyright 2018 C. Leyder


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%INPUT PARAMETERS
% atyp=0;   %analysis type: 0 for 2D, 1 for 3D
% filename  %name of the data file
% 
% flg       %defines the OSP method to be used,
%           %so far available are: efi, efi-dpr, mke and iei (all
%           %sequential algorithms)          
% tm        %Target mode shapes (set to true)
% %ns       %Number of Sensors
% 
% lambda    %Correlation factor %lambda =0.001 produces the uncorrelated model
% ploton=1; %or 0 (=0 => no plots =1 => with plots)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUTS
%Results;      Structure with entries for each final number of sensors and the following fields
%   .pos       Numbers of the DOFs, where a sensor is positioned
%   .Q         Fisher Information matrix of the configuration
%   .PHI       Mode shape vector reduced to entries from the Configuration
%   .H         Information entropy for the configuration
%   .InfoVec   matrix with the information entropy at every iteration step
%   .SensorDof index number of the DOFs to be monitored
%   .sigma     correlation matrix (reduced to the monitored DOFs)
%   .accuracy  determinant of the FIM (reduced to the monitored DOFs)
%   .un         mode shape vector (complete)
%   .EFI        matrix with the effective independence values (for all
%               remaining DOFs)
%   .Qref       %Fischer Information matrix of the reference configuration
%               (full configuration)

%Qref       %Fischer Information matrix of the reference configuration
%           (full configuration)


%% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%MAIN PART
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%Load Mode Shape and Frequency Data from .mat-File
load(filename,'Coord','Modal') 

%% element selection (reduce the mode shape to the selected elements) - not
%necessary for el_sel=0
%el_sel=0 (all), el_sel=1 (columns), el_sel=2 (beams)
if el_sel==1
   a=(Coord.Z==0);  %logical => index whith elements that should be deleted
   %manually keep corner points
   a(2)=0; a(4)=0; a(6)=0; a(8)=0;
   Coord.Z(a)=[];
   Coord.X(a)=[];
   Coord.Y(a)=[];
   Coord.Name(a)=[];
   modeNo = length(Modal.CircFreq);
   for i=1:1:length(a)
   afull(1+modeNo*(i-1):modeNo*i,1)=a(i,1);
   end
   Modal.U1(afull)=[];
   Modal.U3(afull)=[];
   Modal.R2(afull)=[];
elseif el_sel==2
   a=(Coord.Z~=0);  %logical => index which elements should be deleted
   Coord.Z(a)=[];
   Coord.X(a)=[];
   Coord.Y(a)=[];
   Coord.Name(a)=[];
   modeNo = length(Modal.CircFreq);
   for i=1:1:length(a)
   afull(1+modeNo*(i-1):modeNo*i,1)=a(i,1);
   end
   Modal.U1(afull)=[];
   Modal.U3(afull)=[];
   Modal.R2(afull)=[];
end



%% Arrange modal analysis results into a proper mode shape vector (and scale)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%Remark:   
%MM contains 3 columns, one column per DOF (x,y and rot), in each
%column the first 12 displacements are from the first node (one number for
%each mode), the next 12 are for the second node (one number for each mode)
%etc.
%un contains 12 columns, one for each mode
%in each column the first 3 numbers correspond to one node, each number
%corresponds to 1 DOF (x,y, and rot), the next 3 numbers correspond to the
%second node (again one number per DOF (x,y and to)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if atyp==1
MM = [Modal.U1 Modal.U2 Modal.U3 Modal.R1 Modal.R2 Modal.R3];
elseif atyp==0
MM = [Modal.U1 Modal.U3 Modal.R2];  %consider all DOFs%
end

[~,n] = size(MM); 

modeNo = length(Modal.CircFreq);        % number of mode shapes calculated by SAP2000
jointNo = length(Modal.U1)/modeNo;      % nodes of interest defined in SAP, numbered from 1 to 57
dofNo = n*jointNo;                      % number of DOF
un = zeros(modeNo,dofNo);               % vector of mode shapes (transposed)

for ind = 1:jointNo
    un(1:modeNo,(ind-1)*n+1:ind*n) = MM((ind-1)*modeNo+1:ind*modeNo,:);
end

un = un.';
wn = Modal.CircFreq';

for i=1:1:modeNo  %norm all mode shapes to 1 (every mode shape is considered with equal importance) 
    [mval, ind]=max(abs(un(:,i)));
    un(:,i)=un(:,i)*sign(un(ind,i))/mval;
end

% %Calculate the norm of translational and rotational DOFs, so far no
% %special scaling is implemented => add if tiltmeters should be used!
% %split the mode shape vector un into its differnet parts (2 rotational DOFs
% %x and y and one rotational DOF)

%Pre-Allocate
xun=zeros(jointNo,modeNo);
zun=zeros(jointNo,modeNo);
rot2un=zeros(jointNo,modeNo);
normx=zeros(1,modeNo);
normz=zeros(1,modeNo);
normrot=zeros(1,modeNo);

for i=1:1:modeNo 
        xun(1:jointNo,i)=un(1:n:length(un),i);
        zun(1:jointNo,i)=un(2:n:length(un),i);
        rot2un(1:jointNo,i)=un(3:n:(length(un)),i);
        normx(1,i)=norm(xun(:,i));
        normz(1,i)=norm(zun(:,i));
        normrot(1,i)=norm(rot2un(:,i));
        
        un(3:n:(length(un)),i)=un(3:n:(length(un)),i)./1000;
end

%reduce the mode shape vector to the targe mode shapes:
un=un(:,tm);
wn=wn(:,tm);

%%%visualize the mode shapes that are considered in the analysis
if ploton==1
plot_modeshapes(Coord,un,wn,atyp,n)
end

%% Define the correlation error model
d=zeros(dofNo,dofNo); %preallocate

%calculate distances
for i=1:1:dofNo
    for j=1:1:dofNo
d(i,j)=sqrt((Coord.X(ceil(i/3))-Coord.X(ceil(j/3)))^2+(Coord.Z(ceil(i/3))-Coord.Z(ceil(j/3)))^2);
    end
end
%initialize sigma
sigma=eye(dofNo);
%x-Schleife
for i=1:1:dofNo/3
    for j=1:1:dofNo/3
sigma(((i-1)*3+1),((j-1)*3+1))=sqrt(sigma(((i-1)*3+1),((i-1)*3+1))*sigma(((j-1)*3+1),((j-1)*3+1)))*exp(-d(((i-1)*3+1),((j-1)*3+1))/lambda);
    end
end
%z -Schleife
for i=1:1:dofNo/3
    for j=1:1:dofNo/3
sigma(((i-1)*3+2),((j-1)*3+2))=sqrt(sigma(((i-1)*3+2),((i-1)*3+2))*sigma(((j-1)*3+2),((j-1)*3+2)))*exp(-d(((i-1)*3+2),((j-1)*3+2))/lambda);
    end
end
%rot -Schleife
for i=1:1:dofNo/3
    for j=1:1:dofNo/3
sigma(((i-1)*3+3),((j-1)*3+3))=sqrt(sigma(((i-1)*3+3),((i-1)*3+3))*sigma(((j-1)*3+3),((j-1)*3+3)))*exp(-d(((i-1)*3+3),((j-1)*3+3))/lambda);
    end
end

% Sig_maxoffdiag=max(max(sigma-eye(size(sigma))));

Qref=un'*(sigma\un); %Reference FISHER information matrix

%% Perform the OSP-Algorithm
[Results.pos,Results.Q,Results.PHI,Results.H,Results.InfoVec,Results.SensorDof,Results.sigma,Results.accuracy] = OSP_COR_PE_MultiAxis(flg,un,wn,ns,sigma);

%Quality assessment
%information entropy index
Results.Qref=Qref;

%% Plot the information vector (for all DOFs)
% if ploton==1
% figure()
% Color=jet(dofNo);
% hold on;
% 
% % subplot(2,2,1)
% grid on;
% hold on;
% title(strcat('InformationVector-',flg))
% xlabel('DOFs')
% xlim([0 dofNo])
% % upylim=max(max(Results.InfoVec(:,:)));
% % lowylim=min(min(Results.InfoVec(:,:)));
% % ylim([lowylim upylim])
% for i=2+ns:(dofNo-ns+1)
% plot(Results.SensorDof(1:(dofNo-i+1),i-1),Results.InfoVec(1:(dofNo-i+1),i-1),'-*','Color',Color(i,:))
% pause(0.05)
% end
% set(gcf,'paperunits','centimeters')
% set(gcf,'papersize',Papersize) % Desired outer dimensions
% set(gcf,'paperposition',[0,0,Papersize]) % Place plot on figure  
% pause(2)
% 
% % Plot the information vector (on the structure) for x,y and rot DOFs
% % plot InfoVec works only for 2D right now 
% figure()
% %x-direction
% dir=1;
% [Results]=plot_InfoVec(Results,dir,Coord,atyp,ns,length(un(1,:)),jointNo,dofNo,flg);
% set(gcf,'paperunits','centimeters')
% set(gcf,'papersize',Papersize) % Desired outer dimensions
% set(gcf,'paperposition',[0,0,Papersize]) % Place plot on figure  
% 
% 
% 
% figure()
% %y-direction
% dir=2;
% [Results]=plot_InfoVec(Results,dir,Coord,atyp,ns,length(un(1,:)),jointNo,dofNo,flg);
% set(gcf,'paperunits','centimeters')
% set(gcf,'papersize',Papersize) % Desired outer dimensions
% set(gcf,'paperposition',[0,0,Papersize]) % Place plot on figure  
% 
% figure()
% %rot-direction
% dir=3;
% [Results]=plot_InfoVec(Results,dir,Coord,atyp,ns,length(un(1,:)),jointNo,dofNo,flg);
% 
% set(gcf,'paperunits','centimeters')
% set(gcf,'papersize',Papersize) % Desired outer dimensions
% set(gcf,'paperposition',[0,0,Papersize]) % Place plot on figure  
% %print -painters -dpdf -r600 'Z:\01 ETH House of Natural Resources\Presentations\2015 PhD Colloquium Munich\Figures\InfoVec.pdf'

% end
dir=0; %reset to 0 after plotting
%% Plot the final position of sensors (on the structure)

if ploton==1
if atyp==0
    as=0.5; %arrowsize
    al=1.0; %arrowlength   
elseif atyp==1
    as=0.7; %arrowsize
    al=1.0; %arrowlength
end

plot_SensorPos_Multi(Results,Coord,atyp,ns,length(un(1,:)),jointNo,as,al,flg,lambda)
pause(2)
end


% for debugging only
assignin('base', 'un',un);
assignin('base', 'flg',flg);
assignin('base', 'wn',wn);
assignin('base', 'sensorNo',ns);
assignin('base', 'sigma',sigma);
assignin('base', 'd',d);
assignin('base', 'Results',Results);
assignin('base', 'dir',dir);
assignin('base', 'Coord',Coord);
assignin('base', 'jointNo',jointNo);
assignin('base', 'dofNo',dofNo);
end