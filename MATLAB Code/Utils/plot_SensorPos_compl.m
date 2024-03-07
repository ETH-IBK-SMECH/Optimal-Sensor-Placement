function plot_SensorPos_compl(Results,Coord,atyp,ns,~,jointNo,as,al,flg,lambda)
%Function to visualize the final position of the sensors (uni-axial calculation, but complemented with tri-axial sensors in each sensor location)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: C. Leyder
% Last Update: 21.11.2018
% ETH Zurich
% Copyright 2018 C. Leyder 



Fontsize=8;
Width_fig=7.5;
Height_fig=3.;

figure()
hold on;
xlabel('Distance x[m]','Fontsize',Fontsize,'Interpreter','latex');
ylabel('Distance y[m]','Fontsize',Fontsize,'Interpreter','latex');
set(gca,'FontSize',Fontsize,'TickLabelInterpreter','latex');
grid on;

if atyp==0
    n=3;
elseif atyp==1
    n=6;
end

axis equal
    if atyp==1
        view(3)
        set(gca,'XTick',-13:6.5:13)
        set(gca,'YTick',-13:6.5:13)
        set(gca,'ZTick',-4:2:1);
     end

    %set limits for graph
    xmax=max(Coord.X)+1.0;
    xmin=min(Coord.X)-1.0;
    xlim([xmin xmax]);
    
    if atyp==0
    zmax=max(Coord.Z)+1.0;
    zmin=min(Coord.Z)-1.0;
    ylim([zmin zmax]);
    end
    if atyp==1
    xmax=max(Coord.X)+1.0;
    xmin=min(Coord.X)-1.0;
    xlim([xmin xmax]);
    
    ymax=max(Coord.Y)+1.0;
    ymin=min(Coord.Y)-1.0;
    ylim([ymin ymax]);
    
    zmax=max(Coord.Z)+1.0;
    zmin=min(Coord.Z)-1.0;
    zlim([zmin zmax]);   
    end
    
    %axis tight
    box on
   
% titlestr1=sprintf('Sensor Positions (N=%.0f',ns); % 
% titlestr2=sprintf(' & TM=%.0f)',tm);
% titlestr = strcat(titlestr1,titlestr2,'-',flg);
% title(titlestr,'Fontsize',Fontsize+1,'Interpreter','latex')

if strcmp(flg,'efi')
title('Monitored DOFs Positions (EFI Method)','Fontsize',Fontsize+1,'Interpreter','latex');
elseif strcmp(flg,'mke')
title('Monitored DOFs Positions (MKE Method)','Fontsize',Fontsize+1,'Interpreter','latex');
elseif strcmp(flg,'iei') && lambda==0.001
title('Monitored DOFs Positions (IEI Method, $\lambda=0.001$)','Fontsize',Fontsize+1,'Interpreter','latex');
elseif strcmp(flg,'iei') && lambda==0.5
title('Monitored DOFs Positions (IEI Method, $\lambda=0.5$)','Fontsize',Fontsize+1,'Interpreter','latex');
elseif strcmp(flg,'iei') && lambda==1.0
title('Monitored DOFs Positions (IEI Method, $\lambda=1.0$)','Fontsize',Fontsize+1,'Interpreter','latex');
% title('Monitored DOFs Positions (Multi-aixal IEI Method)','Fontsize',Fontsize+1,'Interpreter','latex');
elseif strcmp(flg,'iei') && lambda==1.5
% title('Monitored DOFs Positions (IEI Method, $\lambda=1.5$)','Fontsize',Fontsize+1,'Interpreter','latex');
title('Monitored DOFs (biaxial (complemented))','Fontsize',Fontsize+1,'Interpreter','latex');
end

%UNDEFORMED FRAME
len=length(Coord.Conn(:,1));
X1=zeros(len,1);
X2=zeros(len,1);
Z1=zeros(len,1);
Z2=zeros(len,1);
if atyp==1
Y1=zeros(len,1);
Y2=zeros(len,1);
end

%plot the undeformed frame 
for m=1:1:len
for i=1:1:jointNo
    if strcmp(Coord.Conn{m,1},Coord.Name{i}) 
        X1(m,1)=Coord.X(i);
        if atyp==1
            Y1(m,1)=Coord.Y(i);
        end
        Z1(m,1)=Coord.Z(i);
    end
    if strcmp(Coord.Conn{m,2},Coord.Name{i})
        X2(m,1)=Coord.X(i);
        if atyp==1
            Y2(m,1)=Coord.Y(i);
        end
        Z2(m,1)=Coord.Z(i);
    end
        
end
    if atyp==0
    line([X1(m,1),X2(m,1)],[Z1(m,1),Z2(m,1)],'Color','k');
    elseif atyp==1
    line([X1(m,1),X2(m,1)],[Y1(m,1),Y2(m,1)],[Z1(m,1),Z2(m,1)],'Color','k');
    end
end

    
for i=1:1:ns
    for j=1:1:length(Coord.Z)
        if Results.poscomp(i)==1+n*(j-1)
            if atyp==1
            arrow3d([Coord.X(j) Coord.Y(j) Coord.Z(j)],[Coord.X(j)+al Coord.Y(j) Coord.Z(j)],20,'line',0.3,12)  
            else
            drawarrow_green(Coord.X(j),Coord.X(j)+al,Coord.Z(j),Coord.Z(j),as) 
            end
        elseif Results.poscomp(i)==2+n*(j-1)
            if atyp==1
            arrow3d([Coord.X(j) Coord.Y(j) Coord.Z(j)],[Coord.X(j) Coord.Y(j)+al  Coord.Z(j)],20,'line',0.3,12) 
            else  
            drawarrow_green(Coord.X(j),Coord.X(j),Coord.Z(j),Coord.Z(j)+al,as)   %vertical arrow
            end
        elseif Results.poscomp(i)==3+n*(j-1)
            if atyp==1
            arrow3d([Coord.X(j) Coord.Y(j) Coord.Z(j)],[Coord.X(j) Coord.Y(j)  Coord.Z(j)+al],20,'line',0.3,12) 
            else 
%             line([Coord.X(j),Coord.X(j)],[Coord.Z(j),Coord.Z(j)],'Color',grey,'Marker','o')   %rotational
            end
        elseif atyp==1
            if Results.poscompl(i)==4+n*(j-1)
%               line([Coord.X(j),Coord.X(j)],[Coord.Y(j),Coord.Y(j)],[Coord.Z(j),Coord.Z(j)],'Color',grey,'Marker','o')   
            elseif Results.poscompl(i)==5+n*(j-1)
%               line([Coord.X(j),Coord.X(j)],[Coord.Y(j),Coord.Y(j)],[Coord.Z(j),Coord.Z(j)],'Color',grey,'Marker','o')  
            elseif Results.poscompl(i)==6+n*(j-1)
%               line([Coord.X(j),Coord.X(j)],[Coord.Y(j),Coord.Y(j)],[Coord.Z(j),Coord.Z(j)],'Color',grey,'Marker','o')  
            end
        end           
   end
end

set(gcf,'paperunits','centimeters')
set(gcf,'papersize',[Width_fig,Height_fig])% Desired outer dimensions
set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig])  % Place plot on figure 
if ns==6 && strcmp(flg,'efi')
% print -painters -depsc -r600 'Z:\01 ETH House of Natural Resources\Publications\IALCCE2016\Paper\Figures\SensorPos_EFI.eps'
elseif ns==6 && strcmp(flg,'mke')
% print -painters -depsc -r600 'Z:\01 ETH House of Natural Resources\Publications\IALCCE2016\Paper\Figures\SensorPos_MKE.eps'
elseif (ns==6 && strcmp(flg,'iei') && lambda==0.001)
% print -painters -depsc -r600 'Z:\01 ETH House of Natural Resources\Publications\IALCCE2016\Paper\Figures\SensorPos_IEI_00.eps'
elseif (ns==6 && strcmp(flg,'iei') && lambda==0.5)
% print -painters -depsc -r600 'Z:\01 ETH House of Natural Resources\Publications\IALCCE2016\Paper\Figures\SensorPos_IEI_05.eps'
elseif (ns==6 && strcmp(flg,'iei') && lambda==1.0)
% print -painters -depsc -r600 'Z:\01 ETH House of Natural Resources\Publications\IALCCE2016\Paper\Figures\SensorPos_IEI_10.eps'
elseif (ns==6 && strcmp(flg,'iei') && lambda==1.5)
% print -painters -depsc -r600 'Z:\01 ETH House of Natural Resources\Publications\IALCCE2016\Paper\Figures\SensorPos_IEI_15.eps'
elseif (ns==10 && strcmp(flg,'iei') && lambda==1.5)
% print -painters -depsc -r600 'Z:\01 ETH House of Natural Resources\Publications\IALCCE2016\Paper\Figures\SensorPos_IEI_Uniax.eps'
elseif (ns==21 && strcmp(flg,'iei') && lambda==1.5)
print -painters -depsc -r600 'Z:\01 ETH House of Natural Resources\Publications\IALCCE2016\Paper\Figures\SensorPos_IEI_Uniax_Comp.eps'
end
end