%% PLOTS THE PARETOFRONT FOR DIFFERENT OSP METRICS
% Implemented for dispersion metric (accel, vel and disp), MAC and IEI
% metric

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: C. Leyder
% Last Update: 10.12.2018
% ETH Zurich
% Copyright 2018 C. Leyder

clearvars
close all
clc
grey=[0.65 0.65 0.65];
Fontsize=7;
Width_fig=15;
Height_fig=16;
blueETH=[0.122 0.251 0.478];
Linewidth=1.0;
% Default options for plotting
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'defaultTextFontName', 'Times New Roman')
set(0,'defaultAxesFontSize',Fontsize);
set(0,'defaultTextFontSize',Fontsize);
set(0,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})

load MetricsDispAcc Metrics %load dispersion metrics for acceleration sensors (generated with dispersion/runpadOSP.m)

MetricsN=Metrics;

load MetricsDispVel Metrics  %load dispersion metrics for velocity sensors(generated with dispersion/runpadOSP.m)

MetricsN.DispVel=Metrics.DispVel;

load MetricsDispDisp Metrics %load dispersion metrics for displacement sensors (generated with dispersion/runpadOSP.m)

MetricsN.Dispdisp=Metrics.DispDisp;


load MetricsMAC Metrics     %load MAC metrics (generated with Run_MAC_calculator.m)

MetricsN.MAC=Metrics.MAC;

load IEI_Metric IEI_EFI IEI_IEI IEI_MKE IEI_rand     %load IEI metrics (generated with Comparison_OSP.m and IEI_Metric_Random.m)

MetricsN.IEI.MKE=IEI_MKE;
MetricsN.IEI.EFI=IEI_EFI;
MetricsN.IEI.IEI=IEI_IEI;
MetricsN.IEI.random=IEI_rand;

[~,nsamples]=size(IEI_rand); %number of random samples
nsmax=length(IEI_MKE);
nsmin=find(IEI_MKE>0,1);

sensorRange=nsmin:nsmax;

%% Plot MAC-values for the different modes

figure('DefaultAxesPosition', [0.05, 0.1, 0.9, 0.8])
subplot(6,1,1)
hl1=plot(nsmin:nsmax,MetricsN.MAC.MKE(nsmin:nsmax,1),'k-*');
hold on; grid on;
plot(nsmin:nsmax,MetricsN.MAC.EFI(nsmin:nsmax,1),'-*','Color',grey)
plot(nsmin:nsmax,MetricsN.MAC.IEI(nsmin:nsmax,1),'k-o')
xticks(6:2:28)
xlim([6,28])
ylim([0.7 1.0])
yticks([0.7 0.8 0.9 1.0])
title('Mode 1')
subplot(6,1,2)
plot(nsmin:nsmax,MetricsN.MAC.MKE(nsmin:nsmax,2),'k-*')
hold on; grid on;
hl2=plot(nsmin:nsmax,MetricsN.MAC.EFI(nsmin:nsmax,2),'-*','Color',grey);
plot(nsmin:nsmax,MetricsN.MAC.IEI(nsmin:nsmax,2),'k-o')
title('Mode 2')
xticks(6:2:28)
xlim([6,28])
ylim([0.7 1.0])
yticks([0.7 0.8 0.9 1.0])
subplot(6,1,3)
plot(nsmin:nsmax,MetricsN.MAC.MKE(nsmin:nsmax,3),'k-*')
hold on; grid on;
plot(nsmin:nsmax,MetricsN.MAC.EFI(nsmin:nsmax,3),'-*','Color',grey)
hl3=plot(nsmin:nsmax,MetricsN.MAC.IEI(nsmin:nsmax,3),'k-o');
title('Mode 3')
xticks(6:2:28)
xlim([6,28])
ylim([0.7 1.0])
yticks([0.7 0.8 0.9 1.0])
ylabel('MAC [-]','FontSize',Fontsize);
subplot(6,1,4)
plot(nsmin:nsmax,MetricsN.MAC.MKE(nsmin:nsmax,4),'k-*')
hold on; grid on;
plot(nsmin:nsmax,MetricsN.MAC.EFI(nsmin:nsmax,4),'-*','Color',grey)
plot(nsmin:nsmax,MetricsN.MAC.IEI(nsmin:nsmax,4),'k-o')
title('Mode 4')
xticks(6:2:28)
xlim([6,28])
ylim([0.7 1.0])
yticks([0.7 0.8 0.9 1.0])
subplot(6,1,5)
plot(nsmin:nsmax,MetricsN.MAC.MKE(nsmin:nsmax,5),'k-*');
hold on; grid on;
plot(nsmin:nsmax,MetricsN.MAC.EFI(nsmin:nsmax,5),'-*','Color',grey);
plot(nsmin:nsmax,MetricsN.MAC.IEI(nsmin:nsmax,5),'k-o');
title('Mode 5')
xticks(6:2:28)
xlim([6,28])
ylim([0.7 1.0])
yticks([0.7 0.8 0.9 1.0])
xlabel('\it N_{\rmo}','Fontsize',Fontsize,'Interpreter','tex')
subplot(6,1,6)
axis off
h_legend=legend([hl1,hl2,hl3],'MKE method','EFI method','IEI method');
set(h_legend,'FontSize',Fontsize,'Position',[0.45,0.1,0.1,0.1]);
set(gcf,'papersize',[Width_fig,Height_fig]) % Desired outer dimensions
set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig]) % Place plot on 
% print -painters -dpdf -r600 'MACMetric.pdf'

%% Plot dispersion values, differentiate between acceleration, velocity and displacement sensors
Width_fig=18.6;
Height_fig=7.5;



figure
subplot(1,3,1)
plot(sensorRange,MetricsN.DispAcc.MKE(sensorRange),'k*-') 
axis tight;
hold on; grid on;
plot(sensorRange,MetricsN.DispAcc.EFI(sensorRange),'*-','Color',grey)
plot(sensorRange,MetricsN.DispAcc.IEI(sensorRange),'ko-')
xlabel('Number of monitored DOFs')
ylabel('Energy ratio [%]')
title('Acceleration sensors')
legend('MKE','EFI','IEI','Location','NorthWest')
subplot(1,3,2)
plot(sensorRange,MetricsN.DispVel.MKE(sensorRange),'k*-') 
axis tight;
hold on; grid on;
plot(sensorRange,MetricsN.DispVel.EFI(sensorRange),'*-','Color',grey)
plot(sensorRange,MetricsN.DispVel.IEI(sensorRange),'ko-')
xlabel('Number of monitored DOFs')
ylabel('Energy ratio [%]')
title('Velocity sensors')
% legend('MKE','EFI','IEI')
subplot(1,3,3)
plot(sensorRange,MetricsN.Dispdisp.MKE(sensorRange),'k*-') 
axis tight;
hold on; grid on;
plot(sensorRange,MetricsN.Dispdisp.EFI(sensorRange),'*-','Color',grey)
plot(sensorRange,MetricsN.Dispdisp.IEI(sensorRange),'ko-')
xlabel('Number of monitored DOFs')
ylabel('Energy ratio [%]')
title('Displacement sensors')
% legend('MKE','EFI','IEI')
set(gcf,'papersize',[Width_fig,Height_fig]) % Desired outer dimensions
set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig]) % Place plot on 
% print -painters -dtiff -r600 'DispMetric_Methods.tiff'


figure
subplot(1,3,1)
plot(sensorRange,MetricsN.DispAcc.MKE(sensorRange),'k*-') 
axis tight;
hold on; grid on;
plot(sensorRange,MetricsN.DispVel.MKE(sensorRange),'*-','Color',grey)
plot(sensorRange,MetricsN.Dispdisp.MKE(sensorRange),'ko-')
xlabel('sensor number')
ylabel('Energy ratio (%)')
legend('Acceleration','Velocity','Dispalcement')
title('MKE method')
subplot(1,3,2)
plot(sensorRange,MetricsN.DispAcc.EFI(sensorRange),'k*-') 
axis tight;
hold on; grid on;
plot(sensorRange,MetricsN.DispVel.EFI(sensorRange),'*-','Color',grey)
plot(sensorRange,MetricsN.Dispdisp.EFI(sensorRange),'ko-')
xlabel('sensor number')
ylabel('Energy ratio (%) for velocity sensors')
legend('Acceleration','Velocity','Dispalcement')
title('EFI method')
subplot(1,3,3)
plot(sensorRange,MetricsN.DispAcc.IEI(sensorRange),'k*-') 
axis tight;
hold on; grid on;
plot(sensorRange,MetricsN.DispVel.IEI(sensorRange),'*-','Color',grey)
plot(sensorRange,MetricsN.Dispdisp.IEI(sensorRange),'ko-')
xlabel('sensor number')
ylabel('Energy ratio (%) for velocity sensors')
legend('Acceleration','Velocity','Dispalcement')
title('IEI method')
set(gcf,'papersize',[Width_fig,Height_fig]) % Desired outer dimensions
set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig]) % Place plot on 
% print -painters -dtiff -r600 'Disp_AVD_Metric.tiff'


%% Plot Pareto plots between MAC and dispersion metrics

Width_fig=15;
Height_fig=6;


%build the rms value for the MAC (combine all modes into 1 value)
for i=nsmin:nsmax
MetricsN.MACRMS.MKE(i,1)=rms(MetricsN.MAC.MKE(i,:));
MetricsN.MACRMS.EFI(i,1)=rms(MetricsN.MAC.EFI(i,:));
MetricsN.MACRMS.IEI(i,1)=rms(MetricsN.MAC.IEI(i,:));
    for j=1:1:nsamples
       MetricsN.MACRMS.rand(i,j)=rms(MetricsN.MAC.random(i,j,:));
    end
end

nsl=nsmin;
nsu=nsmax;


%1. Acceleration Sensors (dispersion vs. MAC metric)
XAcc=-[MetricsN.MACRMS.MKE(nsl:nsu),MetricsN.DispAcc.MKE(nsl:nsu)';
      MetricsN.MACRMS.EFI(nsl:nsu),MetricsN.DispAcc.EFI(nsl:nsu)';
      MetricsN.MACRMS.IEI(nsl:nsu),MetricsN.DispAcc.IEI(nsl:nsu)'];
  for i=nsl:1:nsu
        XAcc(end+1:end+nsamples,:)=-[MetricsN.MACRMS.rand(i,:)',MetricsN.DispAcc.rand(i,:)']; 
  end

front = paretofront(XAcc);

figure();
plot(MetricsN.MACRMS.rand(nsl:nsu,:),MetricsN.DispAcc.rand(nsl:nsu,:),'.','Color',grey);
hold on; grid on;
plot(MetricsN.MACRMS.MKE(nsl:nsu),MetricsN.DispAcc.MKE(nsl:nsu),'k*');
plot(MetricsN.MACRMS.EFI(nsl:nsu),MetricsN.DispAcc.EFI(nsl:nsu),'ks');
plot(MetricsN.MACRMS.IEI(nsl:nsu),MetricsN.DispAcc.IEI(nsl:nsu),'ko');
plot(-XAcc(front,1),-XAcc(front,2),'r.','MarkerSize',15);
ylabel('Energy ratio [%]')
xlabel('MAC value [-]')
title('Acceleration sensors')
set(gcf,'papersize',[Width_fig,Height_fig]) % Desired outer dimensions
set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig]) % Place plot on 
% print -painters -dtiff -r600 'Pareto_Acc.tiff'




%2. Velocity sensors
XVel=-[MetricsN.MACRMS.MKE(nsl:nsu),MetricsN.DispVel.MKE(nsl:nsu)';
      MetricsN.MACRMS.EFI(nsl:nsu),MetricsN.DispVel.EFI(nsl:nsu)';
      MetricsN.MACRMS.IEI(nsl:nsu),MetricsN.DispVel.IEI(nsl:nsu)'];
  for i=nsl:nsu
        XVel(end+1:end+nsamples,:)=-[MetricsN.MACRMS.rand(i,:)',MetricsN.DispVel.rand(i,:)']; 
  end

front = paretofront(XVel);

figure();
plot(MetricsN.MACRMS.rand(nsl:nsu,:),MetricsN.DispVel.rand(nsl:nsu,:),'.','Color',grey)
hold on; grid on;
plot(MetricsN.MACRMS.MKE(nsl:nsu),MetricsN.DispVel.MKE(nsl:nsu),'k*')
plot(MetricsN.MACRMS.EFI(nsl:nsu),MetricsN.DispVel.EFI(nsl:nsu),'ks')
plot(MetricsN.MACRMS.IEI(nsl:nsu),MetricsN.DispVel.IEI(nsl:nsu),'ko')
plot(-XVel(front,1),-XVel(front,2),'r.','MarkerSize',15)
ylabel('Energy ratio [%]')
xlabel('MAC value [-]')
title('Velocity sensors')
set(gcf,'papersize',[Width_fig,Height_fig]) % Desired outer dimensions
set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig]) % Place plot on 
% print -painters -dtiff -r600 'Pareto_Vel.tiff'



%3. Displacement sensors

XDisp=-[MetricsN.MACRMS.MKE(nsl:nsu),MetricsN.Dispdisp.MKE(nsl:nsu)';
      MetricsN.MACRMS.EFI(nsl:nsu),MetricsN.Dispdisp.EFI(nsl:nsu)';
      MetricsN.MACRMS.IEI(nsl:nsu),MetricsN.Dispdisp.IEI(nsl:nsu)'];
  for i=nsl:nsu
        XDisp(end+1:end+nsamples,:)=-[MetricsN.MACRMS.rand(i,:)',MetricsN.Dispdisp.rand(i,:)']; 
  end

front = paretofront(XDisp);
figure();
hline1=plot(MetricsN.MACRMS.rand(nsl:nsu,:),MetricsN.Dispdisp.rand(nsl:nsu,:),'.','Color',grey);
hold on; grid on;
hline2=plot(MetricsN.MACRMS.MKE(nsl:nsu),MetricsN.Dispdisp.MKE(nsl:nsu),'k*');
hline3=plot(MetricsN.MACRMS.EFI(nsl:nsu),MetricsN.Dispdisp.EFI(nsl:nsu),'ks');
hline4=plot(MetricsN.MACRMS.IEI(nsl:nsu),MetricsN.Dispdisp.IEI(nsl:nsu),'ko');
hline5=plot(-XDisp(front,1),-XDisp(front,2),'r.','MarkerSize',15);
ylabel('Energy ratio [%]')
xlabel('MAC value [-]')
title('Displacement sensors')
legend([hline1(1,1),hline2(1,1),hline3(1,1),hline4(1,1),hline5(1,1)],'Random configurations','MKE','EFI','IEI','Pareto front','Location','EastOutside')
set(gcf,'papersize',[Width_fig*1.6,Height_fig]) % Desired outer dimensions
set(gcf,'paperposition',[0,0.05,Width_fig*1.6,Height_fig]) % Place plot on 
% print -painters -dtiff -r600 'Pareto_Disp.tiff'



%% Plot Pareto plots between MAC and IEI metric

clear front
XNEW=-[MetricsN.MACRMS.MKE(nsl:nsu),-MetricsN.IEI.MKE(nsl:nsu)';
      MetricsN.MACRMS.EFI(nsl:nsu),-MetricsN.IEI.EFI(nsl:nsu)';
      MetricsN.MACRMS.IEI(nsl:nsu),-MetricsN.IEI.IEI(nsl:nsu)'];
  for i=nsl:1:nsu
        XNEW(end+1:end+nsamples,:)=-[MetricsN.MACRMS.rand(i,:)',-MetricsN.IEI.random(i,:)']; 
  end

front = paretofront(XNEW);

[parplot1,indsort]=sort(-XNEW(front,1));
parplot2=XNEW(front,2);
parplot2=parplot2(indsort);


%Full Pareto Plot
figure();
hline1=semilogy(MetricsN.MACRMS.rand(nsl:nsu,:),MetricsN.IEI.random(nsl:nsu,:),'.','Color',grey);
hold on; grid on;
hline2=semilogy(MetricsN.MACRMS.MKE(nsl:nsu),MetricsN.IEI.MKE(nsl:nsu),'*','Color',[0.4,0.4,0.4]);
hline3=semilogy(MetricsN.MACRMS.EFI(nsl:nsu),MetricsN.IEI.EFI(nsl:nsu),'s','Color',[0.4,0.4,0.4]);
hline4=semilogy(MetricsN.MACRMS.IEI(nsl:nsu),MetricsN.IEI.IEI(nsl:nsu),'o','Color',[0.4,0.4,0.4]);
hline5=semilogy(parplot1,parplot2,'r-','LineWidth',1.5);
legend([hline1(1,1),hline2(1,1),hline3(1,1),hline4(1,1),hline5(1,1)],'Random configurations','MKE','EFI','IEI','Pareto front','Location','EastOutside')
% ylim([-1 10])
ylabel('IEI [-]','Fontsize',Fontsize)
xlabel('RMS(MAC)[-]','Fontsize',Fontsize)
% title('Pareto front','Fontsize',Fontsize)
set(gca,'FontSize',Fontsize);
set(gcf,'papersize',[Width_fig,Height_fig]) % Desired outer dimensions
set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig]) % Place plot on 
% print -painters -dpdf -r600 'OSP_Pareto_full.pdf'


%Zoom-in on the pareto front
figure();
hline1=semilogy(MetricsN.MACRMS.rand(nsl:nsu,:),MetricsN.IEI.random(nsl:nsu,:),'.','Color',grey);
hold on; grid on;
hline2=semilogy(MetricsN.MACRMS.MKE(nsl:nsu),MetricsN.IEI.MKE(nsl:nsu),'*','Color',[0.4,0.4,0.4]);
hline3=semilogy(MetricsN.MACRMS.EFI(nsl:nsu),MetricsN.IEI.EFI(nsl:nsu),'s','Color',[0.4,0.4,0.4]);
hline4=semilogy(MetricsN.MACRMS.IEI(nsl:nsu),MetricsN.IEI.IEI(nsl:nsu),'o','Color',[0.4,0.4,0.4]);
hline5=semilogy(parplot1,parplot2,'r-','LineWidth',0.8);
legend([hline1(1,1),hline2(1,1),hline3(1,1),hline4(1,1),hline5(1,1)],'Random configurations','MKE','EFI','IEI','Pareto front','Location','EastOutside')
ylim([0 10])
ylabel('IEI [-]','Fontsize',Fontsize)
xlabel('RMS(MAC)[-]','Fontsize',Fontsize)
% title('Pareto front','Fontsize',Fontsize)
set(gca,'FontSize',Fontsize);
set(gcf,'papersize',[Width_fig,Height_fig]) % Desired outer dimensions
set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig]) % Place plot on 
% print -painters -dpdf -r600 'OSP_Pareto_zoom.pdf'

%Addition of the implemented configuration (22 DOFs)
figure();
hline1=semilogy(MetricsN.MACRMS.rand(22,:),MetricsN.IEI.random(22,:),'.','Color',grey);
hold on; grid on;
hline2=semilogy(MetricsN.MACRMS.MKE(22),MetricsN.IEI.MKE(22),'*','Color',[0.4,0.4,0.4]);
hline3=semilogy(MetricsN.MACRMS.EFI(22),MetricsN.IEI.EFI(22),'s','Color',[0.4,0.4,0.4]);
hline4=semilogy(MetricsN.MACRMS.IEI(22),MetricsN.IEI.IEI(22),'o','Color',[0.4,0.4,0.4]);
hline5=semilogy(0.996006454105658,2.92309640857134,'k.','MarkerSize',15);
% hline5=semilogy(-XNEW(front,1),XNEW(front,2),'r.','MarkerSize',15);
legend([hline1(1,1),hline2(1,1),hline3(1,1),hline4(1,1),hline5(1,1)],'Random configurations','MKE','EFI','IEI','Implemented configuration','Location','EastOutside')
ylabel('IEI [-]','Fontsize',Fontsize)
ylim([10^0 10^8])
xlim([0.3 1.0])
xlabel('RMS(MAC)[-]','Fontsize',Fontsize)
title('\it N_{\rmo}\rm{=22}','Fontsize',Fontsize,'Interpreter','tex')
set(gca,'FontSize',Fontsize);
set(gcf,'papersize',[Width_fig,Height_fig]) % Desired outer dimensions
set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig]) % Place plot on 
% print -painters -dpdf -r600 'OSP_Config_22.pdf'


%Addition of the implemented configuration (14 DOFs)
figure();
hline1=semilogy(MetricsN.MACRMS.rand(14,:),MetricsN.IEI.random(14,:),'.','Color',grey);
hold on; grid on;
hline2=semilogy(MetricsN.MACRMS.MKE(14),MetricsN.IEI.MKE(14),'*','Color',[0.4,0.4,0.4]);
hline3=semilogy(MetricsN.MACRMS.EFI(14),MetricsN.IEI.EFI(14),'s','Color',[0.4,0.4,0.4]);
hline4=semilogy(MetricsN.MACRMS.IEI(14),MetricsN.IEI.IEI(14),'o','Color',[0.4,0.4,0.4]);
hline5=semilogy(0.996006454105658,3.03032474027788,'k.','MarkerSize',15);
% hline5=semilogy(-XNEW(front,1),XNEW(front,2),'r.','MarkerSize',15);
legend([hline1(1,1),hline2(1,1),hline3(1,1),hline4(1,1),hline5(1,1)],'Random configurations','MKE','EFI','IEI','Implemented configuration','Location','EastOutside')
ylabel('IEI [-]','Fontsize',Fontsize)
ylim([10^0 10^8])
xlim([0.3 1.0])
xlabel('RMS(MAC)[-]','Fontsize',Fontsize)
title('\it N_{\rmo}\rm{=14}','Fontsize',Fontsize,'Interpreter','tex')
set(gca,'FontSize',Fontsize);
set(gcf,'papersize',[Width_fig,Height_fig]) % Desired outer dimensions
set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig]) % Place plot on 
% print -painters -dpdf -r600 'OSP_Config_14.pdf'




