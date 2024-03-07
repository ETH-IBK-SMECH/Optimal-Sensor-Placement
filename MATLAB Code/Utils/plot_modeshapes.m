function  [] = plot_modeshapes(Coord,un,wn,atyp,n)
%%%Function to plot the mode shape from the FEM model (stored in un) 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: C. Leyder
% Last Update: 21.11.2018
% ETH Zurich
% Copyright 2018 C. Leyder


%%% Function INPUTS:
%Coord: Structure with the following fields:
% .X (Vector of X-Coordinates)
% .Y (Vector of Y-Coordinates)
% .Z (Vector of Z-Coordinates)
%. Name (Name of the points)
% .Conn (Matrix with Names of connected points e.g. [A1 A2
%                                                    A2 A3] etc.

%un: Mode Shape Vector (2D structure (atyp=0), [un1x un1y un2x un2y etc.])
%                       (3D structure (atyp=1), [un1x un1y un1z un2x un2y un2z etc.]

%wn: Frequencies corresponding to the mode shape vector
%atyp: Analysis type (0: 2D, 1:3D)
%n: Number of DOFs per point (2: 2D, 3: 3D)


%Graphic definitions
Fontsize=7;
Width_fig=7.5;
Height_fig=3.;
grey=[0.65 0.65 0.65];

% Default options for plotting
set(0,'defaulttextinterpreter','none');
set(0,'defaultAxesFontName', 'Times New Roman')
set(0,'defaultTextFontName', 'Times New Roman')
set(0,'defaultAxesFontSize',Fontsize);
set(0,'defaultTextFontSize',Fontsize);
set(0,{'DefaultAxesXColor','DefaultAxesYColor','DefaultAxesZColor'},{'k','k','k'})


jointNo=length(Coord.X);




%prepare the mode shapes for plotting (separate the different directions)
%Pre-Allocate
xun=zeros(jointNo,length(un(1,:)));
yun=zeros(jointNo,length(un(1,:)));
zun=zeros(jointNo,length(un(1,:)));

for i=1:1:length(un(1,:))
    if atyp==0
            xun(1:jointNo,i)=un(1:n:length(un),i);
            zun(1:jointNo,i)=un(2:n:length(un),i);
    elseif atyp==1
            xun(1:jointNo,i)=un(1:n:length(un),i);
            yun(1:jointNo,i)=un(2:n:length(un),i);
            zun(1:jointNo,i)=un(3:n:length(un),i);
    end
end


for k=1:1:length(un(1,:))
    
    figure()
    hold on;
    set(gca,'FontSize',Fontsize);
    xlabel('x[m]','Fontsize',Fontsize);
    ylabel('y[m]','Fontsize',Fontsize);
    if atyp==1
       zlabel('Distance z[m]','Fontsize',Fontsize);
    end
    
    freq=wn(k)/(2*pi());
    fvalue=sprintf('=%.1f',freq);
    titel=strcat({'Mode Shape '},num2str(k),{',\it{  f} \rm\bf'},fvalue,{' [Hz]'});
    title(titel,'Fontsize',Fontsize,'Interpreter','tex');
    grid on;
    axis equal;
    box on;
    
    if atyp==1
    view(3)
    end

    if atyp==1
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
    elseif atyp==1
    ymax=max(Coord.Y)+1.5;
    ymin=min(Coord.Y)-1.5;
    ylim([ymin ymax]);
    zmax=max(Coord.Z)+1.5;
    zmin=min(Coord.Z)-1.5;
    zlim([zmin zmax]);
    end
    
%scaling factor for the plotting of mode shapes
if atyp==1
sc=20;
elseif atyp==0
sc=1;
end

for i=1:1:jointNo
Coorddef.X(i,k)=Coord.X(i)+sc*xun(i,k);
    if atyp==1
    Coorddef.Y(i,k)=Coord.Y(i)+sc*yun(i,k);
    end
Coorddef.Z(i,k)=Coord.Z(i)+sc*zun(i,k);
end

%predefine vectors
len=length(Coord.Conn(:,1));
X1=zeros(len,1);
X2=zeros(len,1);
Z1=zeros(len,1);
Z2=zeros(len,1);
if atyp==1
Y1=zeros(len,1);
Y2=zeros(len,1);
end


%% plot the unddeformed frame
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

%% Plot the deformed frame
%predefine vectors
X1def=zeros(len,1);
X2def=zeros(len,1);
if atyp==1
Y1def=zeros(len,1);
Y2def=zeros(len,1);
end
Z1def=zeros(len,1);
Z2def=zeros(len,1);

%plot the mode shapes 
for m=1:1:len
for i=1:1:jointNo
    if strcmp(Coord.Conn{m,1},Coord.Name{i}) 
        X1def(m,1)=Coorddef.X(i,k);
        if atyp==1
           Y1def(m,1)=Coorddef.Y(i,k); 
        end
        Z1def(m,1)=Coorddef.Z(i,k);
    end
    if strcmp(Coord.Conn{m,2},Coord.Name{i})
        X2def(m,1)=Coorddef.X(i,k);
        if atyp==1
            Y2def(m,1)=Coorddef.Y(i,k);
        end
        Z2def(m,1)=Coorddef.Z(i,k);
    end
        
end
if atyp==0
line([X1def(m,1),X2def(m,1)],[Z1def(m,1),Z2def(m,1)],'Color',grey); 
elseif atyp==1
    line([X1def(m,1),X2def(m,1)],[Y1def(m,1),Y2def(m,1)],[Z1def(m,1),Z2def(m,1)],'Color',grey); 
end
end

set(gcf,'paperunits','centimeters')
set(gcf,'papersize',[Width_fig,Height_fig])% Desired outer dimensions
set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig])  % Place plot on figure 
if k==1
print -painters -dpdf -r600 'ModeShape1.pdf'
elseif k==2
print -painters -dpdf -r600 'ModeShape2.pdf'
elseif k==3
print -painters -dpdf -r600 'ModeShape3.pdf'   
elseif k==4
print -painters -dpdf -r600 'ModeShape4.pdf'
elseif k==5
print -painters -dpdf -r600 'ModeShape5.pdf'  
elseif k==6
print -painters -dpdf -r600 'ModeShape6.pdf' 
elseif k==7
print -painters -dpdf -r600 'ModeShape7.pdf' 
elseif k==8
print -painters -dpdf -r600 'ModeShape8.pdf' 
end

end