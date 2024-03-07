function  [MAC] = plot_reduced_modeshapes(Coord,un,wn,atyp,n,sens_ind,ploton)

%%%Function to interpolate and plot the mode shapes. Additionally the
%%%MAC-value between the original mode shape from the FEM model and an 
%interpolated mode shape (based on the reduced information available from 
%the sensor configuration) (=the MAC metric) is evaluated.

%NOTE: For the interpolation several assumptions are necessary (e.g.
%columns and beams do not elongate, etc.). These assumptions need to be
%adjusted for each problem setting.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: C. Leyder
% Last Update: 5.12.2018
% ETH Zurich
% Copyright 2018 C. Leyder


Fontsize=8;
grey=[0.65 0.65 0.65];

jointNo=length(Coord.X);


% Order the reference FE DOFs
    if atyp==0
            xun = un(1:n:end,:);
            zun = un(2:n:end,:);
    elseif atyp==1
            xun = un(1:n:end,:);
            yun = un(2:n:end,:);
            zun = un(3:n:end,:);
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
empty_ind.x = setdiff(1:jointNo,sens_ind.x);    %this is the index of non-sensor positions along x
empty_ind.z = setdiff(1:jointNo,sens_ind.z);    %this is the index of non-sensor positions along z

% Group the nodes per span V: vertical, H:horizontal, # from left to right
V1 = find(Coord.X==-9.75);  %column 1
V2 = find(Coord.X==-3.25);  %column 2   
V3 = find(Coord.X==3.25);   %column 3
V4 = find(Coord.X==9.75);   %column 4

H1 = find(Coord.X>=-9.75&Coord.X<=-3.25&Coord.Z==0);    %beam 1
H2 = find(Coord.X>=-3.25&Coord.X<=3.25&Coord.Z==0);     %beam 2
H3 = find(Coord.X>=3.25&Coord.X<=9.75&Coord.Z==0);      %beam 3



Coordint.X = repmat(Coord.X,size(un,2)); %Initialize coordinates
Coordint.Z = repmat(Coord.Z,size(un,2)); %Initialize coordinates

for k=1:1:length(un(1,:))
    
if ploton==1 
    figure()
    hold on;
    set(gca,'FontSize',Fontsize,'TickLabelInterpreter','latex');
    xlabel('Distance x[m]','Fontsize',Fontsize,'Interpreter','latex');
    ylabel('Distance y[m]','Fontsize',Fontsize,'Interpreter','latex');
    if atyp==1
       zlabel('Distance z[m]','Fontsize',Fontsize);
    end
    
    freq=wn(k)/(2*pi());
    fvalue=sprintf('=%.1f',freq);
    titel=strcat({'Mode Shape '},num2str(k),{', f'},fvalue,{' [Hz]'});
    title(titel,'Fontsize',Fontsize+1,'Interpreter','latex');
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
    
end
% scaling factor for the plotting of mode shapes
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
if ploton==1 

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
%     line([X1(m,1),X2(m,1)],[Z1(m,1),Z2(m,1)],'Color','k');
    elseif atyp==1
%     line([X1(m,1),X2(m,1)],[Y1(m,1),Y2(m,1)],[Z1(m,1),Z2(m,1)],'Color','k');
    end
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

%% Plot the mode shapes 
if ploton==1 

for m=1:len
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
line([X1(m,1),X1def(m,1)],[Z1(m,1),Z1def(m,1)],'Color',grey); 
line([X2(m,1),X2def(m,1)],[Z2(m,1),Z2def(m,1)],'Color',grey); 
elseif atyp==1
line([X1def(m,1),X2def(m,1)],[Y1def(m,1),Y2def(m,1)],[Z1def(m,1),Z2def(m,1)],'Color',grey); 
end
end

end

%% Interpolating at Positions different to the Sensors
temps=[];
tempe=[];
% Span V1
for i=1
sV1 = intersect(V1,[sens_ind.x; 1]);   %Location 1 has 0 displacement in x & y
if i==1
Coordint.X(sV1,k) = Coorddef.X(sV1,k);  %Coordint.X was initialized with the undeformed coordinates
% if isempty(intersect(sV1,2))
%     Coordint.X(2,k) = Coord.X(2)+mean(Coorddef.X(intersect(H1,[sens_ind.x]))-Coord.X(intersect(H1,[sens_ind.x])));  %assume that the top node moves along with the beam
%     sV1(end+1)=2;
% end

else
Coordint.Z(sV1,k) = Coorddef.Z(sV1,k);
end

eV1 = intersect(V1,empty_ind.x);
eV1 = setdiff(eV1,1);    

temps = [temps;sV1];
tempe = [tempe;eV1];
if length(sV1)>1
    if i==1
Coordint.X(eV1,k) = interp1(Coord.Z(sV1),Coorddef.X(sV1,k),Coord.Z(eV1),'pchip');
    else
Coordint.Z(eV1,k) = interp1(Coord.Z(sV1),Coorddef.Z(sV1,k),Coord.Z(eV1),'pchip');
    end
end
end

if ploton==1 
scatter(Coorddef.X(temps,k),Coorddef.Z(temps,k),'b')
scatter(Coordint.X(tempe,k),Coordint.Z(tempe,k),'r');
end

temps=[];
tempe=[];
% Span V2
for i=1
sV2 = intersect(V2,[sens_ind.x; 3]);   %Location 3 has 0 displacement in x & y; Location 3 corresponds to DOFind=5
if i==1
Coordint.X(sV2,k) = Coorddef.X(sV2,k);
% % % % if isempty(intersect(sV2,4))
% % % %     Coordint.X(4,k) = Coord.X(4)+mean([mean(Coorddef.X(intersect(H1,[sens_ind.x]))-Coord.X(intersect(H1,[sens_ind.x]))),mean(Coorddef.X(intersect(H2,[sens_ind.z]))-Coord.X(intersect(H2,[sens_ind.z])))]);  %assume that the top node moves along with the beam average left and right
% % % %     sV2(end+1)=4;
% % % % end
else
Coordint.Z(sV2,k) = Coorddef.Z(sV2,k);
end

eV2 = intersect(V2,empty_ind.x);
eV2 = setdiff(eV2,3);

temps = [temps;sV2];
tempe = [tempe;eV2];
if length(sV2)>1
    if i==1
Coordint.X(eV2,k) = interp1(Coord.Z(sV2),Coorddef.X(sV2,k),Coord.Z(eV2),'pchip');
    else
Coordint.Z(eV2,k) = interp1(Coord.Z(sV2),Coorddef.Z(sV2,k),Coord.Z(eV2),'pchip');
    end
end
end

if ploton==1 
scatter(Coorddef.X(temps,k),Coorddef.Z(temps,k),'b')
scatter(Coordint.X(tempe,k),Coordint.Z(tempe,k),'r');
end

temps=[];
tempe=[];
% Span V3
for i=1
sV3 = intersect(V3,[sens_ind.x; 5]);   %Location 5 has 0 displacement in x & y (corresponds to DOFind = 7)
if i==1
Coordint.X(sV3,k) = Coorddef.X(sV3,k);

% if isempty(intersect(sV3,6))
%     Coordint.X(6,k) = Coord.X(6)+mean([mean(Coorddef.X(intersect(H2,[sens_ind.x]))-Coord.X(intersect(H2,[sens_ind.x]))),mean(Coorddef.X(intersect(H3,[sens_ind.z]))-Coord.X(intersect(H3,[sens_ind.z])))]);  %assume that the top node moves along with the beam average left and right
%     sV3(end+1)=6;
% end


else
Coordint.Z(sV3,k) = Coorddef.Z(sV3,k);
end

eV3 = intersect(V3,empty_ind.x);
eV3 = setdiff(eV3,5);

temps = [temps;sV3];
tempe = [tempe;eV3];
if length(sV3)>1
    if i==1
Coordint.X(eV3,k) = interp1(Coord.Z(sV3),Coorddef.X(sV3,k),Coord.Z(eV3),'pchip');
    else
Coordint.Z(eV3,k) = interp1(Coord.Z(sV3),Coorddef.Z(sV3,k),Coord.Z(eV3),'pchip');
    end
end
end
% 
if ploton==1 
scatter(Coorddef.X(temps,k),Coorddef.Z(temps,k),'b')
scatter(Coordint.X(tempe,k),Coordint.Z(tempe,k),'r');
end

temps=[];
tempe=[];
% Span V4
for i=1
sV4 = intersect(V4,[sens_ind.x; 7]);   %Location 7 has 0 displacement in x & y (Coordind = 3)
if i==1
Coordint.X(sV4,k) = Coorddef.X(sV4,k);

% if isempty(intersect(sV4,8))
%     Coordint.X(8,k) = Coord.X(8)+mean(Coorddef.X(intersect(H3,[sens_ind.x]))-Coord.X(intersect(H3,[sens_ind.x])));  %assume that the top node moves along with the beam average left and right
%     sV4(end+1)=8;
% end

else
Coordint.Z(sV4,k) = Coorddef.Z(sV4,k);
end

eV4 = intersect(V4,empty_ind.x);
eV4 = setdiff(eV4,7);

temps = [temps;sV4];
tempe = [tempe;eV4];
if length(sV4)>1
    if i==1
Coordint.X(eV4,k) = interp1(Coord.Z(sV4),Coorddef.X(sV4,k),Coord.Z(eV4),'pchip');
    else
Coordint.Z(eV4,k) = interp1(Coord.Z(sV4),Coorddef.Z(sV4,k),Coord.Z(eV4),'pchip');
    end
end
end

if ploton==1 
scatter(Coorddef.X(temps,k),Coorddef.Z(temps,k),'b')
scatter(Coordint.X(tempe,k),Coordint.Z(tempe,k),'r');
end
    
temps=[];
tempe=[];
% Span H1
l1=2;
l2=4;
for i=2
sH1 = intersect(H1,[sens_ind.z; l1;l2]);   %Locations 2,4 have 0 displacement in y 
if i==1
Coordint.X(sH1,k) = Coorddef.X(sH1,k);
else
Coordint.Z(sH1,k) = Coorddef.Z(sH1,k);
Coordint.Z(l1,k) = Coord.Z(l1);   %Locations 2,4 have 0 displacement in y
Coordint.Z(l2,k) = Coord.Z(l2);
sH1_x = intersect(H1,[sens_ind.x]);
Coordint.X(sH1_x,k) = Coorddef.X(sH1_x,k);
end

%eH1 contains the nodes with the unmeasured coordinates
eH1 = intersect(H1,empty_ind.z);
eH1 = setdiff(eH1,[l1 l2]);

    
temps = [temps;sH1];
tempe = [tempe;eH1];
if length(sH1)>1
    if i==1
Coordint.X(eH1,k) = interp1(Coord.X(sH1),Coorddef.X(sH1,k),Coord.X(eH1),'spline');
    else
        if isempty(intersect(H1,[sens_ind.x]))==false
          Coordint.X([eH1;intersect(H1,sens_ind.z)],k) = Coord.X([eH1;intersect(H1,sens_ind.z)])+interp1(Coord.X(intersect(H1,[sens_ind.x;l1;l2])),Coordint.X(intersect(H1,[sens_ind.x;l1;l2]),k)-Coord.X(intersect(H1,[sens_ind.x;l1;l2])),Coord.X([eH1;intersect(H1,sens_ind.z)]),'linear');
           else
          Coordint.X([eH1;intersect(H1,sens_ind.z)],k) = Coord.X([eH1;intersect(H1,sens_ind.z)])+interp1(Coord.X([l1,l2]),Coordint.X([l1,l2],k)-Coord.X([l1,l2]),Coord.X([eH1;intersect(H1,sens_ind.z)]),'linear');  
        end
Coordint.Z(eH1,k) = interp1(Coord.X(sH1),Coorddef.Z(sH1,k),Coord.X(eH1),'spline');
    end
end
end

if ploton==1 
scatter(Coordint.X(temps,k),Coordint.Z(temps,k),'b')  %measured / known points
scatter(Coordint.X([l1;l2],k),Coordint.Z([l1;l2],k),'r') 
scatter(Coordint.X(tempe,k),Coordint.Z(tempe,k),'r');
scatter(Coordint.X(sH1_x,k),Coordint.Z(sH1_x,k),'b');
end

temps=[];
tempe=[];
% Span H2
l1=4;
l2=6;
for i=2
sH2 = intersect(H2,[sens_ind.z; l1;l2]);   %Locations 4,6 have 0 displacement in y 
if i==1
Coordint.X(sH2,k) = Coorddef.X(sH2,k);
else
Coordint.Z(sH2,k) = Coorddef.Z(sH2,k);
Coordint.Z(l1,k) = Coord.Z(l1);   %Locations 4,6 have 0 displacement in y
Coordint.Z(l2,k) = Coord.Z(l2);  %add x sensors knowledge
sH2_x = intersect(H2,[sens_ind.x]);
Coordint.X(sH2_x,k) = Coorddef.X(sH2_x,k);
end

eH2 = intersect(H2,empty_ind.z);
eH2 = setdiff(eH2,[l1 l2]);

temps = [temps;sH2];
tempe = [tempe;eH2];
if length(sH2)>1
    if i==1
Coordint.X(eH2,k) = interp1(Coord.X(sH2),Coorddef.X(sH2,k),Coord.X(eH2),'pchip');
    else
        if isempty(intersect(H2,[sens_ind.x]))==false
         Coordint.X([eH2;intersect(H2,sens_ind.z)],k) = Coord.X([eH2;intersect(H2,sens_ind.z)])+interp1(Coord.X(intersect(H2,[sens_ind.x;l1;l2])),Coordint.X(intersect(H2,[sens_ind.x;l1;l2]),k)-Coord.X(intersect(H2,[sens_ind.x;l1;l2])),Coord.X([eH2;intersect(H2,sens_ind.z)]),'linear');      %average known x deformations (assume beam =rigid body motion in x direction)
        else
          Coordint.X([eH2;intersect(H2,sens_ind.z)],k) = Coord.X([eH2;intersect(H2,sens_ind.z)])+interp1(Coord.X([l1,l2]),Coordint.X([l1,l2],k)-Coord.X([l1,l2]),Coord.X([eH2;intersect(H2,sens_ind.z)]),'linear');  
        end
    Coordint.Z(eH2,k) = interp1(Coord.X(sH2),Coorddef.Z(sH2,k),Coord.X(eH2),'pchip');
    end
end
end

if ploton==1 
scatter(Coordint.X(temps,k),Coordint.Z(temps,k),'b')
scatter(Coordint.X([l1;l2],k),Coordint.Z([l1;l2],k),'r') 
scatter(Coordint.X(tempe,k),Coordint.Z(tempe,k),'r');
scatter(Coordint.X(sH2_x,k),Coordint.Z(sH2_x,k),'b');
end

temps=[];
tempe=[];
% Span H3
l1=6;
l2=8;
for i=2
sH3 = intersect(H3,[sens_ind.z; l1;l2]);   %Locations 6,8 have 0 displacement in y 
if i==1
Coordint.X(sH3,k) = Coorddef.X(sH3,k);
else
Coordint.Z(sH3,k) = Coorddef.Z(sH3,k);
Coordint.Z(l1,k) = Coord.Z(l1);   %Locations 6,8 have 0 displacement in y
Coordint.Z(l2,k) = Coord.Z(l2);
sH3_x = intersect(H3,[sens_ind.x]);
Coordint.X(sH3_x,k) = Coorddef.X(sH3_x,k);
end

%eH3 contains the nodes with the unmeasured coordinates
eH3 = intersect(H3,empty_ind.z);
eH3 = setdiff(eH3,[l1 l2]);

temps = [temps;sH3];
tempe = [tempe;eH3];
if length(sH3)>1
    if i==1
Coordint.X(eH3,k) = interp1(Coord.X(sH3),Coorddef.X(sH3,k),Coord.X(eH3),'pchip');
    else
        if isempty(intersect(H3,[sens_ind.x]))==false
          Coordint.X([eH3;intersect(H3,sens_ind.z)],k) = Coord.X([eH3;intersect(H3,sens_ind.z)])+interp1(Coord.X(intersect(H3,[sens_ind.x;l1;l2])),Coordint.X(intersect(H3,[sens_ind.x;l1;l2]),k)-Coord.X(intersect(H3,[sens_ind.x;l1;l2])),Coord.X([eH3;intersect(H3,sens_ind.z)]),'linear');      %average known x deformations (assume beam =rigid body motion in x direction)
        else
          Coordint.X([eH3;intersect(H3,sens_ind.z)],k) = Coord.X([eH3;intersect(H3,sens_ind.z)])+interp1(Coord.X([l1,l2]),Coordint.X([l1,l2],k)-Coord.X([l1,l2]),Coord.X([eH3;intersect(H3,sens_ind.z)]),'linear');  
        end
    Coordint.Z(eH3,k) = interp1(Coord.X(sH3),Coorddef.Z(sH3,k),Coord.X(eH3),'pchip');
    end
end
end

if ploton==1 
scatter(Coordint.X(temps,k),Coordint.Z(temps,k),'b')
scatter(Coordint.X([l1;l2],k),Coordint.Z([l1;l2],k),'r') 
scatter(Coordint.X(tempe,k),Coordint.Z(tempe,k),'r');
scatter(Coordint.X(sH3_x,k),Coordint.Z(sH3_x,k),'b');
end

% % Plot the Spans in groups
% scatter(Coord.X(V1),Coord.Z(V1),'r')
% scatter(Coord.X(V2),Coord.Z(V2),'b')
% scatter(Coord.X(V3),Coord.Z(V3),'r')
% scatter(Coord.X(V4),Coord.Z(V4),'b')
% scatter(Coord.X(H1),Coord.Z(H1),'g')
% scatter(Coord.X(H2),Coord.Z(H2),'m')
% scatter(Coord.X(H3),Coord.Z(H3),'g')

% set(gcf,'paperunits','centimeters')
% set(gcf,'papersize',[Width_fig,Height_fig])% Desired outer dimensions
% set(gcf,'paperposition',[0,0.05,Width_fig,Height_fig])  % Place plot on figure 
% if k==1
% % print -painters -dpdf -r600 'SAPModeShape1.pdf'
% elseif k==2
% % print -painters -dpdf -r600 'SAPModeShape2.pdf'
% elseif k==3
% % print -painters -dpdf -r600 'SAPModeShape3.pdf'   
% elseif k==4
% % print -painters -dpdf -r600 'SAPModeShape4.pdf'
% elseif k==5
% % print -painters -dpdf -r600 'SAPModeShape5.pdf'  
% elseif k==6
% % print -painters -dpdf -r600 'SAPModeShape6.pdf' 
% elseif k==7
% % print -painters -dpdf -r600 'SAPModeShape7.pdf' 
% elseif k==8
% % print -painters -dpdf -r600 'SAPModeShape8.pdf' 
% end


end

%% Calculate the MAC Values
mode_int=zeros(64*2,size(un,2));
mode_ref=mode_int;
MAC=zeros(1,size(un,2));
for k=1:size(un,2)
mode_int(:,k)=[Coordint.X(1:64,k)-Coord.X(1:64);
    Coordint.Z(1:64,k)-Coord.Z(1:64);];
mode_ref(:,k)=[Coorddef.X(1:64,k)-Coord.X(1:64);
    Coorddef.Z(1:64,k)-Coord.Z(1:64);];
MAC(k)=(abs(mode_int(:,k)'*mode_ref(:,k)))^2/((mode_int(:,k)'*mode_int(:,k))*(mode_ref(:,k)'*mode_ref(:,k)));
end
