%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements optimal sensor placement algorithms for multi-axial sensors, 
% with considering their multi-axiality during the optimization process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: C. Leyder
% Last Update: 20.11.2018
% ETH Zurich
% Copyright 2018 C. Leyder

clearvars;
clc;
close all;
startanalysis=tic;
disp(tic);

Fontsize=8;
Width_fig=7.5;
Height_fig=6;
grey=[0.65 0.65 0.65];

addpath('./Utils')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%INPUT PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
atyp=0;   %analysis type: 0 for 2D, 1 for 3D
filename='DataExcel_2DFrame.mat';

%OSP Algorithm
%flg='efi';  %defines the OSP method to be used,
            %so far available are: efi, efi-dpr, mke and bssp (all
            %sequential algorithms)
%flags={'efi';'efi-dpr';'mke';'bssp'};  %use all 4 models
flags={'iei'};
%Target mode shapes             
tm=false(12,1);   %mode shapes to be considered => put value to 1!
tm(1:5)=1;
%Number of Sensors
% ns=length(nonzeros(tm))+1; 
% ns=9:3:15;
    %vary number of sensors
    nsmin=9;
    nsmax=42;
%Correlation factor
% lambda=0.001;
% lambda=[0.001 0.5 1.0 1.5];  %lambda =0 produces the uncorrelated model, tested also lambda =2.0, the iei is then reduced!
lambda=1.5;
%element selection
el_sel=0; %el_sel=0 (all), el_sel=1 (columns), el_sel=2 (beams)
%%%%%%with plots or without (1=on, 0=off)
ploton=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Function Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. Main File: Comparison_OSP_Multi.m
%2. Functions in Main File: Evaluation_OSP_Multi.m;
%3. Functions in Evaluation_OSP.m:
%                               plot_modeshapes.m
%                               OSP_COR_PE.m
%                               plot_InfoVec.m
%                               plot_SensorPos.m
%4. Functions in plot_SensorPos.m: arrow3d.m and drawarrow_green.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%Execute OSP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Pre-Allocate:
%Metrics for configurations evaluated with the EFI method (Information
%entropy index (IEI), determinant of the Fischer Information Matrix (detQ))
IEI_EFI=zeros(length(lambda),length(nsmax-nsmin));
detQ_EFI=zeros(length(lambda),length(nsmax-nsmin));
%idem for configurations derived with the MKE method
IEI_MKE=zeros(length(lambda),length(nsmax-nsmin));
detQ_MKE=zeros(length(lambda),length(nsmax-nsmin));
%idem for configurations derived with the IEI method
IEI_IEI=zeros(length(lambda),length(nsmax-nsmin));
detQ_IEI=zeros(length(lambda),length(nsmax-nsmin));

detQref=zeros(length(lambda),1); % determinant of the fischer information matrix of the reference configuration

for indm=1:1:length(flags)  %index methods
    method=flags(indm);
    disp(method)
for ns=nsmin:3:nsmax %index 3-axial sensors (step=3)
    numsens = ns;
    disp(numsens)
for lm=1:1:length(lambda)
    lam=lambda(lm);
    disp(lam)
[Ns_Results(ns),Qref]=Evaluation_OSP_Multi(flags{indm},filename,atyp,tm,ns,lambda(lm),el_sel,ploton); %#ok<SAGROW>
    if strcmp(flg,'efi')
    IEI_EFI(lm,ns)=sqrt(det(Qref)/det(Ns_Results(ns).Q));
    detQ_EFI(lm,ns)=det(Ns_Results(ns).Q);
    elseif strcmp(flg,'mke')
    IEI_MKE(lm,ns)=sqrt(det(Qref)/det(Ns_Results(ns).Q));
    detQ_MKE(lm,ns)=det(Ns_Results(ns).Q);
    elseif strcmp(flg,'iei')
    IEI_IEI(lm,ns)=sqrt(det(Qref)/det(Ns_Results(ns).Q));
    detQ_IEI(lm,ns)=det(Ns_Results(ns).Q);
    end
    detQref(lm)=det(Qref);
end
end
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%Compare Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot the prediction error correlation factor versus the determinant of
% the reference configuration (shows that this value is dependent on the
% value of lambda and can therefore not be selected for direct comparisons)
% figure()
% plot(lambda,detQref)
% hold on;
% xlabel('\lambda');
% ylabel('det(Qref)');
% grid on;


% Information enthropy visualisation (compare configurations from MKE and
% IEI)
% 
% figure(); 
% x=nsmin:1:nsmax;
% h1=semilogy(x,IEI_MKE(1,nsmin:1:nsmax),'k-*'); 
% hold on;
% h2=semilogy(x,IEI_EFI(1,nsmin:1:nsmax),'-*','Color',grey);
% grid on;
% xlabel('Number of monitored DOFs','Fontsize',Fontsize,'Interpreter','latex')
% ylabel('IEI [-]','Fontsize',Fontsize,'Interpreter','latex')
% set(gca,'FontSize',Fontsize,'TickLabelInterpreter','latex');
% title('Information entropy index','Fontsize',Fontsize+1,'Interpreter','latex')
% h_legend=legend([h1 h2],{['MKE']; ['EFI']},'location','NorthEast');
% set(h_legend,'FontSize',Fontsize-1,'Interpreter','latex');
% set(gcf,'paperunits','centimeters')
% set(gcf,'papersize',[Width_fig,Height_fig])% Desired outer dimensions
% set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig])  % Place plot on figure 
% print -painters -depsc -r600 'IEI_EFIANDMKE.eps'



%% compare different lambda values (IEI method)

% figure(); 
% x=nsmin:1:nsmax;
% h(1)=semilogy(x,IEI_IEI(1,nsmin:1:nsmax),'k-*');
% hold on;
% grid on;
% h(2)=semilogy(x,IEI_IEI(2,nsmin:1:nsmax),'-*','Color',grey);
% h(3)=semilogy(x,IEI_IEI(3,nsmin:1:nsmax),'k-o');
% h(4)=semilogy(x,IEI_IEI(4,nsmin:1:nsmax),'-o','Color',grey);
% xlabel('Number of monitored DOFs','Fontsize',Fontsize,'Interpreter','latex')
% ylabel('IEI [-]','Fontsize',Fontsize,'Interpreter','latex')
% set(gca,'FontSize',Fontsize,'TickLabelInterpreter','latex');
% % ylim([10^0  10^6])
% title('Information entropy index','Fontsize',Fontsize+1,'Interpreter','latex')
% h_legend=legend([h(1) h(2) h(3) h(4)],{['$\lambda$ = 0.001']; ['$\lambda$ = 0.5']; ['$\lambda$ = 1.0']; ['$\lambda$ = 1.5']},'location','NorthEast');
% set(h_legend,'FontSize',Fontsize-1,'Interpreter','latex');
% set(gcf,'paperunits','centimeters')
% set(gcf,'papersize',[Width_fig,Height_fig])% Desired outer dimensions
% set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig])  % Place plot on figure 
% print -painters -depsc -r600 'Z:\01 ETH House of Natural Resources\Publications\IALCCE2016\Paper\Figures\IEI_Lambdas.eps'

%% compare single axis / multi-axis
%gather data from multi-axis
Multi=nsmin:3:nsmax;
Multi(2,:)=IEI_IEI(nsmin:3:nsmax);
Multi(3,:)=2/3.*Multi(1,:);  %for the 2D case, only 2 of the 3 available DOFs are monitored at each position

figure(); 
h(2)=semilogy(Multi(3,:),Multi(2,:),'k-*');
hold on;
grid on;
xlabel('Number of monitored DOFs','Fontsize',Fontsize,'Interpreter','latex')
ylabel('IEI [-]','Fontsize',Fontsize,'Interpreter','latex')
set(gca,'FontSize',Fontsize,'TickLabelInterpreter','latex');
title('Information entropy index','Fontsize',Fontsize+1,'Interpreter','latex')
h_legend=legend(h(2),{'Multi-axis'},'location','NorthEast');
set(h_legend,'FontSize',Fontsize-1,'Interpreter','latex');
set(gcf,'paperunits','centimeters')
set(gcf,'papersize',[Width_fig,Height_fig])% Desired outer dimensions
set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig])  % Place plot on figure 
% print -painters -depsc -r600 'IEI_Lambdas_biaxial.eps'

save('MultiAxis.mat','Multi');


toc(startanalysis)