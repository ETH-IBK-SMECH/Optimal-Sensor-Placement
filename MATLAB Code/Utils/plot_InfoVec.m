function [Results]=plot_InfoVec(Results,dir,Coord,atyp,ns,tm,jointNo,dofNo,flg)
%Function to plot the evaluation of the information vector across the
%iteration steps for each direction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: C. Leyder
% Last Update: 21.11.2018
% ETH Zurich
% Copyright 2018 C. Leyder

xlabel('Distance x[m]');
ylabel('Distance y[m]');
grid on;
Color=jet(dofNo);


axis equal
    if atyp==1
        view(3)
        set(gca,'XTick',-13:6.5:13)
        set(gca,'YTick',-13:6.5:13)
        set(gca,'ZTick',-4:2:1);
     end

    %set limits for graph
    xmax=max(Coord.X)+2.5;
    xmin=min(Coord.X)-2.5;
    xlim([xmin xmax]);
    
    if atyp==0
    zmax=max(Coord.Z)+5.5;
    zmin=min(Coord.Z)-1.0;
    ylim([zmin zmax]);
    end
    if atyp==1
    xmax=max(Coord.X)+2.5;
    xmin=min(Coord.X)-2.5;
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
   
    if dir==1
        titlestr1=sprintf('x-dir (N=%.0f',ns); 
    elseif dir==2
        titlestr1=sprintf('y-dir (N=%.0f',ns); 
    elseif dir==3
        titlestr1=sprintf('rot-dir (N=%.0f',ns); 
    end
titlestr2 = sprintf(' & TM=%.0f)',tm);
titlestr = strcat(titlestr1,titlestr2,'-',flg);
title(titlestr)

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
    line([X1(m,1),X2(m,1)],[Z1(m,1),Z2(m,1)]);
    elseif atyp==-1 || atyp==-2  
    line([X1(m,1),X2(m,1)],[Z1(m,1),Z2(m,1)]);
    elseif atyp==1
    line([X1(m,1),X2(m,1)],[Y1(m,1),Y2(m,1)],[Z1(m,1),Z2(m,1)]);
    end
end
hold on;
%PLOT THE INFORMATION VECTOR - ON THE FRAME
%choose 1 DOF, either rot, x or y (separately!!)

%create a new matrix, that contains the x and z-coordinates of the DOFs to
%be monitored
for p=1:1:(dofNo-ns)
for j=1:1:dofNo/3
a=find(Results.SensorDof(:,p)==1+j*3-3); %x-DOF
b=find(Results.SensorDof(:,p)==2+j*3-3); %z-DOF
c=find(Results.SensorDof(:,p)==3+j*3-3); %rot-DOF

if dir==1
  %x-dir = value of interest
  Results.Coord.X(a,p)=Coord.X(j,1); 
  Results.Coord.Z(a,p)=Coord.Z(j,1);
  %set the 2 other directions to NaN
  Results.Coord.X(b,p)=NaN; 
  Results.Coord.Z(b,p)=NaN; 
  Results.Coord.X(c,p)=NaN; 
  Results.Coord.Z(c,p)=NaN;
  
  
elseif dir==2
  %x-dir = value of interest
  Results.Coord.X(b,p)=Coord.X(j,1); 
  Results.Coord.Z(b,p)=Coord.Z(j,1);
  %set the 2 other directions to NaN
  Results.Coord.X(a,p)=NaN; 
  Results.Coord.Z(a,p)=NaN; 
  Results.Coord.X(c,p)=NaN; 
  Results.Coord.Z(c,p)=NaN;  

elseif dir==3
  %x-dir = value of interest
  Results.Coord.X(c,p)=Coord.X(j,1); 
  Results.Coord.Z(c,p)=Coord.Z(j,1);
  %set the 2 other directions to NaN
  Results.Coord.X(a,p)=NaN; 
  Results.Coord.Z(a,p)=NaN; 
  Results.Coord.X(b,p)=NaN; 
  Results.Coord.Z(b,p)=NaN;  
        
end
end
   
end

ilb=6;
iub=(dofNo-ns);

%scaling factor

Results.sc=5.0/max(Results.InfoVec(:,end));
Results.sc=1.0;
Results.shift(:,1)=Results.sc*Results.InfoVec(:,1);



 for i=ilb:1:iub
%           plot(Results.Coord.X(1:(dofNo-i+1),i-1),Results.Coord.Z(1:(dofNo-i+1),i-1)+Results.sc*Results.InfoVec(1:(dofNo-i+1),i-1)-Results.shift(1:(dofNo-i+1),1),'*','Color',Color(i,:));
          plot(Results.Coord.X(1:(dofNo-i+1),i-1),Results.Coord.Z(1:(dofNo-i+1),i-1)+Results.sc*Results.InfoVec(1:(dofNo-i+1),i-1),'*','Color',Color(i,:));
          pause(0.05);
 end
 pause(1)

 
end