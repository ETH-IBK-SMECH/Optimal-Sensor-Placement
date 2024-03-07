function [r_ind]=randomConfig(nsensors,nsamples,dofNo)

%Generates random configurations by permuting nsensors from an available
%set of DOFs (dofNo), nsamples-times

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Author: C. Leyder
% Last Update: 4.12.2018
% ETH Zurich
% Copyright 2018 C. Leyder

for i=1:1:nsamples
        r_ind(1:nsensors,i)=randperm(dofNo,nsensors);  %random samples with unique values (sampling without replacement!!!)
end
end