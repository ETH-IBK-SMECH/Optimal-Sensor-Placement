%Generate Random Configurations
% Generates nsamples of random sensor configurations for a set number of
% dof (It is assumed that for a 2D problem, only the x and y DOFs can be
% monitored, the x-DOFs sensors are positioned on DOFs 1,4,7, etc.
% the y-DOFs sensors are positioned on DOFs 2,5,8, etc.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: C. Leyder
% Last Update: 4.12.2018
% ETH Zurich
% Copyright 2018 C. Leyder

clearvars; close all; clc;
dofNo=128;   % from the OSP methods no rotational DOFs are monitored => throw these out
             % total DOF No=192, x,y only! (But DOF numbering is 1(x),2(y),3(rot), etc. 
nsamples=100; %number of random configurations
nsmin=6;  %minimum number of monitored DOFs
nsmax=28; %maximum number of monitored DOFs
i=0;


randomindexes=zeros(nsmax,nsamples,nsmax);

for nsensors=nsmin:1:nsmax
% nsensors=6;
i=i+1;
randomindexes(1:nsensors,1:nsamples,i)=randomConfig(nsensors,nsamples,dofNo); %sampling without replacement
    for jj=1:1:nsamples
       for kk=1:1:nsensors
           if rem(randomindexes(kk,jj,i),2)==1  %uneven number
               randind(kk,jj,i)=(((randomindexes(kk,jj,i)+1)/2)-1)*3+1;  %translate to DOF numbering 
           elseif rem(randomindexes(kk,jj,i),2)==0 %even number
               randind(kk,jj,i)=(((randomindexes(kk,jj,i))/2)-1)*3+2;   %translate to DOF numbering
           end
       end
     end

end

save('RandomConfigs.mat','randind','randomindexes')