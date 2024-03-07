%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implements optimal sensor placement algorithms for uni-axial sensors, or
% multi-axial sensors, without considering their multi-axiality
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
disp(tic)

Fontsize=7;
Width_fig=7.5;
Height_fig=6;
grey=[0.65 0.65 0.65];


% Default options for plotting
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'defaultTextFontName', 'Times New Roman')
set(0,'defaultAxesFontSize',Fontsize);
set(0,'defaultTextFontSize',Fontsize);
set(0,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

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
flags={'efi','mke','iei'};
%Target mode shapes             
tm=false(12,1);   %mode shapes to be considered => put value to 1!
tm(1:5)=1;
%Number of Sensors
    %define the range of the total number of sensors (from nsmin to nsmax)
    nsmin=6;
    nsmax=28;
%Prediction error correlation factor (only affects the IEI method, not the
%EFI and MKE method)!
% lambda=[0.001,0.5,1.0,1.5];  %lambda =0 produces the uncorrelated model, tested also lambda =2.0, the iei is then reduced!
lambda=1.5;
%element selection
el_sel=0; %el_sel=0 (all), el_sel=1 (columns), el_sel=2 (beams)
%%%%%%with plots or without (1=on, 0=off)
ploton=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%Function Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1. Main File: Comparison_OSP.m
%2. Function in Main File: Evaluation_OSP.m;
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

%Pre-allocate:

%Evaluation metrics for configurations derived with the EFI method
%(Information entropy index (IEI), determinant of the FIM (detQ), Maximum
%Kinetic Energy value, Effective independence value, trace of the FIM
%(trace_Q), Eigenvector value product
IEI_EFI=zeros(length(lambda),length(nsmax-nsmin));
detQ_EFI=zeros(length(lambda),length(nsmax-nsmin));
MKE_EFI=zeros(length(lambda),length(nsmax-nsmin));
EFI_EFI=zeros(length(lambda),length(nsmax-nsmin));
Trace_Q_EFI=zeros(length(lambda),length(nsmax-nsmin));
EVP_EFI=zeros(length(lambda),length(nsmax-nsmin));
%Evaluation metrics for configurations derived with the MKE method
IEI_MKE=zeros(length(lambda),length(nsmax-nsmin));
detQ_MKE=zeros(length(lambda),length(nsmax-nsmin));
MKE_MKE=zeros(length(lambda),length(nsmax-nsmin));
EFI_MKE=zeros(length(lambda),length(nsmax-nsmin));
Trace_Q_MKE=zeros(length(lambda),length(nsmax-nsmin));
EVP_MKE=zeros(length(lambda),length(nsmax-nsmin));
%Evaluation metrics for configurations derived with the IEI method
IEI_IEI=zeros(length(lambda),length(nsmax-nsmin));
detQ_IEI=zeros(length(lambda),length(nsmax-nsmin));
MKE_IEI=zeros(length(lambda),length(nsmax-nsmin));
EFI_IEI=zeros(length(lambda),length(nsmax-nsmin));
Trace_Q_IEI=zeros(length(lambda),length(nsmax-nsmin));
EVP_IEI=zeros(length(lambda),length(nsmax-nsmin));

detQref=zeros(length(lambda),1); %determinant of the FIM of the reference onfiguration, the reference configuration is the full configuration, where each DOF is equipped with a sensor
EVP_ii=zeros(length(nsmax-nsmin)); %eigenvector value product for a single configuration with a fixed number of sensors

for indm=1:1:length(flags)  %index methods
    method=flags(indm);
    disp(method)
for ns=nsmin:1:nsmax
    numsens = ns;
    disp(ns)
for lm=1:1:length(lambda)
    lam=lambda(lm);
    disp(lam)
[Ns_Results(ns),Qref]=Evaluation_OSP(flags{indm},filename,atyp,tm,ns,lambda(lm),el_sel,ploton); %#ok<SAGROW>
    if strcmp(flg,'efi')
    IEI_EFI(lm,ns)=sqrt(det(Qref)/det(Ns_Results(ns).Q));
    detQ_EFI(lm,ns)=det(Ns_Results(ns).Q);
    MKE_EFI(lm,ns)=sum(sum(Ns_Results(ns).PHI.*Ns_Results(ns).PHI,2))/sum(sum(Ns_Results(ns).un.*Ns_Results(ns).un,2));
    EFI_EFI(lm,ns)=sum(Ns_Results(ns).InfoVec(1:ns,end))/sum(Ns_Results(ns).InfoVec(:,1));
    Trace_Q_EFI(lm,ns)=trace(Ns_Results(ns).Q);
    for ii=1:1:ns
    EVP_ii(ii)=prod(Ns_Results(ns).PHI(ii,:));
    end
    EVP_EFI(lm,ns)=sum(EVP_ii);
    elseif strcmp(flg,'mke')
    IEI_MKE(lm,ns)=sqrt(det(Qref)/det(Ns_Results(ns).Q));
    detQ_MKE(lm,ns)=det(Ns_Results(ns).Q);
    MKE_MKE(lm,ns)=sum(Ns_Results(ns).InfoVec(1:ns,end))/sum(Ns_Results(ns).InfoVec(:,1));
    EFI_MKE(lm,ns)=sum(Ns_Results(ns).EFI(1:ns,end))/sum(Ns_Results(ns).EFI(:,1));
    Trace_Q_MKE(lm,ns)=trace(Ns_Results(ns).Q);
    for ii=1:1:ns
    EVP_ii(ii)=prod(Ns_Results(ns).PHI(ii,:));
    end
    EVP_MKE(lm,ns)=sum(EVP_ii);
    elseif strcmp(flg,'iei')
    IEI_IEI(lm,ns)=sqrt(det(Qref)/det(Ns_Results(ns).Q));
    detQ_IEI(lm,ns)=det(Ns_Results(ns).Q);
    MKE_IEI(lm,ns)=sum(sum(Ns_Results(ns).PHI.*Ns_Results(ns).PHI,2))/sum(sum(Ns_Results(ns).un.*Ns_Results(ns).un,2));
    EFI_IEI(lm,ns)=sum(Ns_Results(ns).EFI(1:ns,end))/sum(Ns_Results(ns).EFI(:,1));
    Trace_Q_IEI(lm,ns)=trace(Ns_Results(ns).Q);
    for ii=1:1:ns
    EVP_ii(ii)=prod(Ns_Results(ns).PHI(ii,:));
    end
    EVP_IEI(lm,ns)=sum(EVP_ii);
    end
    detQref(lm)=det(Qref);
end
end
if strcmp(flg,'efi')
save('./mac calculation/Ns_Results_EFI.mat','Ns_Results')  %save here for evaluation of the Modal Assurance Criterion Metric
elseif strcmp(flg,'mke')
save('./mac calculation/Ns_Results_MKE.mat','Ns_Results')  %save here for evaluation of the Modal Assurance Criterion Metric
elseif strcmp(flg,'iei')
save('./mac calculation/Ns_Results_IEI.mat','Ns_Results')  %save here for evaluation of the Modal Assurance Criterion Metric
end
end


save('./paretofront/IEI_Metric.mat','IEI_IEI','IEI_MKE','IEI_EFI')  %save here for inclusion of the IEI metric in the ParetoFront calculations

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


% % Information enthropy visualisation (compare MKE and EFI)
% % 
% figure(); 
% x=nsmin:1:nsmax;
% h1=semilogy(x,IEI_MKE(1,nsmin:1:nsmax),'k-*'); 
% hold on;
% h2=semilogy(x,IEI_EFI(1,nsmin:1:nsmax),'-*','Color',grey);
% h3=semilogy(x,IEI_IEI(1,nsmin:1:nsmax),'-*','Color',[0.4 0.4 0.4]);
% grid on;
% xlabel('\it N_{\rmo}','Fontsize',Fontsize,'Interpreter','tex')
% ylabel('IEI [-]','Fontsize',Fontsize)
% set(gca,'FontSize',Fontsize);
% % title('Information entropy index','Fontsize',Fontsize)
% h_legend=legend([h1 h2 h3],{['MKE']; ['EFI']; ['IEI (\lambda=0.001)']},'location','NorthEast');
% set(h_legend,'FontSize',Fontsize);
% set(gcf,'paperunits','centimeters')
% set(gcf,'papersize',[Width_fig,Height_fig])% Desired outer dimensions
% set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig])  % Place plot on figure 
% % print -painters -dpdf -r600 'IEI_EFIANDMKE.pdf'


% %% compare different lambda values (IEI method)
% figure(); 
% x=nsmin:1:nsmax;
% h(1)=semilogy(x,IEI_IEI(1,nsmin:1:nsmax),'k-*');
% hold on;
% grid on;
% h(2)=semilogy(x,IEI_IEI(2,nsmin:1:nsmax),'-*','Color',grey);
% h(3)=semilogy(x,IEI_IEI(3,nsmin:1:nsmax),'k-o');
% h(4)=semilogy(x,IEI_IEI(4,nsmin:1:nsmax),'-o','Color',grey);
% xlabel('\it N_{\rmo}','Fontsize',Fontsize,'Interpreter','tex')
% ylabel('IEI [-]','Fontsize',Fontsize)
% set(gca,'FontSize',Fontsize);
% % ylim([10^0  10^6])
% % title('Information entropy index','Fontsize',Fontsize)
% h_legend=legend([h(1) h(2) h(3) h(4)],{['\lambda_{cor} = 0.001']; ['\lambda_{cor} = 0.5']; ['\lambda_{cor} = 1.0']; ['\lambda_{cor} = 1.5']},'location','NorthEast');
% set(h_legend,'FontSize',Fontsize);
% set(gcf,'paperunits','centimeters')
% set(gcf,'papersize',[Width_fig,Height_fig])% Desired outer dimensions
% set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig])  % Place plot on figure 
% % print -painters -dpdf -r600 'IEI_Lambdas.pdf'


% PLOT MKE vs IEI metrics
% figure()
% %MKE vs. IEI values for sensors obtained via the EFI method
% semilogy(1-MKE_EFI(nsmin:nsmax),IEI_EFI(nsmin:nsmax),'*-')
% hold on; grid on;
% %MKE vs. IEI values for sensors obtained via the MKE method
% semilogy(1-MKE_MKE(nsmin:nsmax),IEI_MKE(nsmin:nsmax),'*-')
% %MKE vs. IEI values for sensors obtained via the IEI method
% semilogy(1-MKE_IEI(nsmin:nsmax),IEI_IEI(nsmin:nsmax),'*-')
% legend('EFI','MKE','IEI','Location','Eastoutside')
% xlabel('Maximum Kinetic Energy')
% ylabel('Information Entropy Index')
% semilogy([1-MKE_IEI(6),1-MKE_EFI(6),1-MKE_MKE(6)],[IEI_IEI(6),IEI_EFI(6),IEI_MKE(6)],'k-')
% if nsmax>9
% semilogy([1-MKE_IEI(10),1-MKE_EFI(10),1-MKE_MKE(10)],[IEI_IEI(10),IEI_EFI(10),IEI_MKE(10)],'k-')
% end
% if nsmax>19
% semilogy([1-MKE_IEI(20),1-MKE_EFI(20),1-MKE_MKE(20)],[IEI_IEI(20),IEI_EFI(20),IEI_MKE(20)],'k-')
% end

% PLOT MKE vs EFI Metric
% figure()
% %MKE vs. EFI values for sensors obtained via the EFI method
% semilogy(1-MKE_EFI(nsmin:nsmax),1-EFI_EFI(nsmin:nsmax),'*-')
% hold on; grid on;
% %MKE vs. EFI values for sensors obtained via the MKE method
% semilogy(1-MKE_MKE(nsmin:nsmax),1-EFI_MKE(nsmin:nsmax),'*-')
% %MKE vs. EFI values for sensors obtained via the IEI method
% semilogy(1-MKE_IEI(nsmin:nsmax),1-EFI_IEI(nsmin:nsmax),'*-')
% legend('EFI','MKE','IEI','Location','Eastoutside')
% xlabel('Maximum Kinetic Energy')
% ylabel('Sum(Effective Independence)')
% semilogy([1-MKE_IEI(6),1-MKE_EFI(6),1-MKE_MKE(6)],[1-EFI_IEI(6),1-EFI_EFI(6),1-EFI_MKE(6)],'k-')
% if nsmax>9
% semilogy([1-MKE_IEI(10),1-MKE_EFI(10),1-MKE_MKE(10)],[1-EFI_IEI(10),1-EFI_EFI(10),1-EFI_MKE(10)],'k-')
% end
% if nsmax>19
% semilogy([1-MKE_IEI(20),1-MKE_EFI(20),1-MKE_MKE(20)],[1-EFI_IEI(20),1-EFI_EFI(20),1-EFI_MKE(20)],'k-')
% end
% figure()
% plot(nsmin:nsmax,Trace_Q_EFI(nsmin:nsmax),'k-*')
% hold on; grid on;
% plot(nsmin:nsmax,Trace_Q_MKE(nsmin:nsmax),'-*','Color',grey)
% plot(nsmin:nsmax,Trace_Q_IEI(nsmin:nsmax),'k-o')
% h=legend('EFI','MKE','IEI','Location','SouthEast');
% title('Comparison of OSP methods','Fontsize',Fontsize+1,'Interpreter','latex')
% xlabel('Number of monitored DOFs','Fontsize',Fontsize,'Interpreter','latex')
% ylabel('Trace of Q [-]','Fontsize',Fontsize,'Interpreter','latex')
% set(gca,'FontSize',Fontsize,'TickLabelInterpreter','latex');
% set(h,'FontSize',Fontsize-1,'Interpreter','latex');
% set(gcf,'paperunits','centimeters')
% set(gcf,'papersize',[Width_fig,Height_fig])% Desired outer dimensions
% set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig])  % Place plot on figure 
% % print -painters -dpdf -r600 'Z:\01 ETH House of Natural Resources\Presentations\2016_10_IALCCE\Figures\TraceQ.pdf'
% 

% Plot EVP (Eigenvector value product)
% figure()
% plot(nsmin:nsmax,EVP_EFI(nsmin:nsmax))
% hold on; grid on;
% plot(nsmin:nsmax,EVP_MKE(nsmin:nsmax))
% plot(nsmin:nsmax,EVP_IEI(nsmin:nsmax))
% legend('EFI','MKE','IEI','Location','Eastoutside')
% title('Eigenvalue vector product')
% xlabel('Number of monitored DOFs','Fontsize',Fontsize,'Interpreter','latex')
% ylabel('Eigenvalue of vector product')

%% Prepare data for comparison single and multi-axis
%%% gather data from single-axis
Single=nsmin:1:nsmax;
Single(2,:)=IEI_IEI(nsmin:1:nsmax);

%complement (for the case that everywhere that a sensor should be
%positioned according to the 1-axial case, in fact a tri-axial sensor is
%placed (here times 2 as we are considering the 2D case)
Single(3,:)=2.*Single(1,:);

%Pre-allocate
IEI_IEI_comp=zeros(length(nsmax-nsmin),1);


for n=nsmin:1:nsmax
    Ns_Results(n).poscomp=zeros(length(Ns_Results(n).pos)*3,1);
    j=1;
for i=1:1:length(Ns_Results(n).pos)
    if j==1 || ( (j>1) && (max(Ns_Results(n).poscomp(1:j-1))<Ns_Results(n).pos(i)))
        if rem(Ns_Results(n).pos(i)-1,3)==0  %x-dir
        Ns_Results(n).poscomp(j)=Ns_Results(n).pos(i);
        Ns_Results(n).poscomp(j+1)=Ns_Results(n).pos(i)+1;
        Ns_Results(n).poscomp(j+2)=Ns_Results(n).pos(i)+2;
        elseif rem(Ns_Results(n).pos(i)-1,3)==1 %y-dir
        Ns_Results(n).poscomp(j)=Ns_Results(n).pos(i)-1;
        Ns_Results(n).poscomp(j+1)=Ns_Results(n).pos(i);
        Ns_Results(n).poscomp(j+2)=Ns_Results(n).pos(i)+1;
        elseif rem(Ns_Results(n).pos(i)-1,3)==2 %rot-dir
        Ns_Results(n).poscomp(j)=Ns_Results(n).pos(i)-2;
        Ns_Results(n).poscomp(j+1)=Ns_Results(n).pos(i)-1;
        Ns_Results(n).poscomp(j+2)=Ns_Results(n).pos(i);   
        end
    end
    j=j+3;
end
%compute the FIM for the complemented setup (multi-axial sensors positioned
%on spots, where a single axis sensor should be positioned according to the
%OSP)
[L,U] = lu(sigma(nonzeros(Ns_Results(n).poscomp),nonzeros(Ns_Results(n).poscomp)));
Ns_Results(n).Qcomp=un(nonzeros(Ns_Results(n).poscomp),:)'*(U\(L\un(nonzeros(Ns_Results(n).poscomp),:)));
Ns_Results(n).Qcheck=un(Ns_Results(n).pos,:)'*(sigma(Ns_Results(n).pos,Ns_Results(n).pos)\un(Ns_Results(n).pos,:));
Ns_Results(n).Check=round(Ns_Results(n).Qcheck\Ns_Results(n).Q,1);

IEI_IEI_comp(n)=sqrt(det(Ns_Results(n).Qref)/det(Ns_Results(n).Qcomp));
end

Single(4,:)=IEI_IEI_comp(nsmin:nsmax);


save('SingleAxis.mat','Single');

toc(startanalysis)

% save('./paretofront/OSP_Results.mat') %save here for further processing