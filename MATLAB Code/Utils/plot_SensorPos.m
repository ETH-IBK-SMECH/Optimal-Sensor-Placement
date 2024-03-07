function plot_SensorPos(Results,Coord,atyp,ns,~,jointNo,as,al,flg,lambda)
%Function to visualize the final position of the sensors (uni-axial)

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
xlabel('x[m]','Fontsize',Fontsize);
ylabel('y[m]','Fontsize',Fontsize);
set(gca,'FontSize',Fontsize);
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
    else
        set(gca,'XTick',-13:6.5:13)
        set(gca,'YTick',-4:2:1)
     end

    %set limits for graph
    xmax=max(Coord.X)+1.5;
    xmin=min(Coord.X)-1.5;
    xlim([xmin xmax]);
    
    if atyp==0
    zmax=max(Coord.Z)+1.5;
    zmin=min(Coord.Z)-1.5;
    ylim([zmin zmax]);
    end
    if atyp==1
    xmax=max(Coord.X)+1.5;
    xmin=min(Coord.X)-1.5;
    xlim([xmin xmax]);
    
    ymax=max(Coord.Y)+1.5;
    ymin=min(Coord.Y)-1.5;
    ylim([ymin ymax]);
    
    zmax=max(Coord.Z)+1.5;
    zmin=min(Coord.Z)-1.5;
    zlim([zmin zmax]);   
    end
    
    %axis tight
    box on
   
% titlestr1=sprintf('Sensor Positions (N=%.0f',ns); % 
% titlestr2=sprintf(' & TM=%.0f)',tm);
% titlestr = strcat(titlestr1,titlestr2,'-',flg);
% title(titlestr,'Fontsize',Fontsize+1,'Interpreter','latex')

if strcmp(flg,'efi')
% title('Sensor configuration (EFI)','Fontsize',Fontsize);
elseif strcmp(flg,'mke')
% title('Sensor configuration (MKE)','Fontsize',Fontsize);
elseif strcmp(flg,'iei') && lambda==0.001
title('\lambda_{cor}=0.001','Fontsize',Fontsize,'Interpreter','tex');
elseif strcmp(flg,'iei') && lambda==0.5
title('\lambda_{cor}=0.5','Fontsize',Fontsize,'Interpreter','tex');
elseif strcmp(flg,'iei') && lambda==1.0
title('\lambda_{cor}=1.0','Fontsize',Fontsize,'Interpreter','tex');
% title('Monitored DOFs Positions (Multi-aixal IEI Method)','Fontsize',Fontsize+1,'Interpreter','latex');
elseif strcmp(flg,'iei') && lambda==1.5
title('\lambda_{cor}=1.5','Fontsize',Fontsize,'Interpreter','tex');
% title('Sensor configuration (uniaxial)','Fontsize',Fontsize+1);
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
        if Results.pos(i)==1+n*(j-1)
            if atyp==1
            arrow3d([Coord.X(j) Coord.Y(j) Coord.Z(j)],[Coord.X(j)+al Coord.Y(j) Coord.Z(j)],20,'line',0.3,12)  
            else
            drawarrow_green(Coord.X(j),Coord.X(j)+al,Coord.Z(j),Coord.Z(j),as) 
            end
        elseif Results.pos(i)==2+n*(j-1)
            if atyp==1
            arrow3d([Coord.X(j) Coord.Y(j) Coord.Z(j)],[Coord.X(j) Coord.Y(j)+al  Coord.Z(j)],20,'line',0.3,12) 
            else  
            drawarrow_green(Coord.X(j),Coord.X(j),Coord.Z(j),Coord.Z(j)+al,as)   %vertical arrow
            end
        elseif Results.pos(i)==3+n*(j-1)
            if atyp==1
            arrow3d([Coord.X(j) Coord.Y(j) Coord.Z(j)],[Coord.X(j) Coord.Y(j)  Coord.Z(j)+al],20,'line',0.3,12) 
            else 
%             line([Coord.X(j),Coord.X(j)],[Coord.Z(j),Coord.Z(j)],'Color',grey,'Marker','o')   %rotational
            end
        elseif atyp==1
            if Results.pos(i)==4+n*(j-1)
%               line([Coord.X(j),Coord.X(j)],[Coord.Y(j),Coord.Y(j)],[Coord.Z(j),Coord.Z(j)],'Color',grey,'Marker','o')   
            elseif Results.pos(i)==5+n*(j-1)
%               line([Coord.X(j),Coord.X(j)],[Coord.Y(j),Coord.Y(j)],[Coord.Z(j),Coord.Z(j)],'Color',grey,'Marker','o')  
            elseif Results.pos(i)==6+n*(j-1)
%               line([Coord.X(j),Coord.X(j)],[Coord.Y(j),Coord.Y(j)],[Coord.Z(j),Coord.Z(j)],'Color',grey,'Marker','o')  
            end
        end           
   end
end

set(gcf,'paperunits','centimeters')
set(gcf,'papersize',[Width_fig,Height_fig])% Desired outer dimensions
set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig])  % Place plot on figure 
if ns==6 && strcmp(flg,'efi')
% print -painters -dpdf -r600 'SensorPos_EFI.pdf'
elseif ns==6 && strcmp(flg,'mke')
% print -painters -dpdf -r600 'SensorPos_MKE.pdf'
elseif (ns==6 && strcmp(flg,'iei') && lambda==0.001)
% print -painters -dpdf -r600 'SensorPos_IEI_00.pdf'
elseif (ns==6 && strcmp(flg,'iei') && lambda==0.5)
% print -painters -dpdf -r600 'SensorPos_IEI_05.pdf'
elseif (ns==6 && strcmp(flg,'iei') && lambda==1.0)
% print -painters -dpdf -r600 'SensorPos_IEI_10.pdf'
elseif (ns==6 && strcmp(flg,'iei') && lambda==1.5)
% print -painters -dpdf -r600 'SensorPos_IEI_15.pdf'
elseif (ns==10 && strcmp(flg,'iei') && lambda==1.5)
% print -painters -depsc -r600 'SensorPos_IEI_Uniax.eps'
elseif (ns==7 && strcmp(flg,'iei') && lambda==1.5)
% print -painters -depsc -r600 'SensorPos_IEI_Uniax_Comp.eps'
end
end