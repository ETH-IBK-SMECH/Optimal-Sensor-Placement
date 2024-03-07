%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare configurations evaluted with Comparison_OSP.m and
% Comparison_OSP_Multi.m to show the effect of considering the
% muli-axiality of sensors during the optimization problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: C. Leyder
% Last Update: 20.11.2018
% ETH Zurich
% Copyright 2018 C. Leyder
clearvars
close all
clc

Fontsize=8;
Width_fig=7.5;
Height_fig=6;
grey=[0.65 0.65 0.65];

addpath('./data')

load('SingleAxis.mat');
load('MultiAxis.mat');

figure(); 
h(1)=semilogy(Single(1,:),Single(2,:),'k-*');
hold on;
grid on;
h(2)=semilogy(Single(3,1:end),Single(4,1:end),'-*','Color',grey);
h(3)=semilogy(Multi(3,1:end),Multi(2,1:end),'k-o');
xlabel('Number of monitored DOFs','Fontsize',Fontsize,'Interpreter','latex')
ylabel('IEI [-]','Fontsize',Fontsize,'Interpreter','latex')
set(gca,'FontSize',Fontsize,'TickLabelInterpreter','latex');
title('Information entropy index','Fontsize',Fontsize+1,'Interpreter','latex')
h_legend=legend([h(1) h(2) h(3)],{'uniaxial';'biaxial (complemented)';'biaxial (calculated)'},'location','NorthEast');
set(h_legend,'FontSize',Fontsize-1,'Interpreter','latex');
set(gcf,'paperunits','centimeters')
set(gcf,'papersize',[Width_fig,Height_fig])% Desired outer dimensions
set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig])  % Place plot on figure 
print -painters -depsc -r600 'IEI_biaxial.eps'
