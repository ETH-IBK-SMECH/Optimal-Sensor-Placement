function [Delta,wn,D,d] = dispersionAnalysis(sys,SigmaF)

%
% dispersionAnalysis
%
% performs dispersion analysis of structural systems.
%
% [Delta,wn] = dispersionAnalysis(sys,SigmaF)
%
% returns the percentage modal dispersions DELTA that correspond 
% to vibration modes with frequencies WN (Hz). SYS is the system 
% in state-space format SigmaF the zero lag covariance matrix of 
% the input.
%
% [Delta,wn,D,d] = dispersionAnalysis(sys,SigmaF)
%
% further returns the modal dispersions D in non percentage form
% and the related eigenvalues of the state transition matrix.
%
% dispersionAnalysis(sys,SigmaF)
%
% prints the results on workspace.
%
% See also modalAnalysis, modalTransferNorm, id2da.
%
% Reference page in Help browser:
% <a href="matlab: web([docroot '/toolbox/mdac/funref/dispersionAnalysis.html'],'-helpbrowser')">doc dispersionAnalysis</a>
%

%
% Author: V. Ntertimanis
% 1st Ed: 10-08-2010
% Last Update: 09-03-2017
% National Technical University of Athens
% Copyright 1995-2013 V.K. Ntertimanis
%

if ~ isa(sys,'ss')
    sys = ss(sys);
end
% state - space matrices
Ass = sys.A; n = size(Ass,1);
Bss = sys.B;
Css = sys.C; m = size(Css,1);
Dss = sys.D;
Ts = sys.Ts;
% required matrices
SIGMA = Bss*SigmaF*Bss';
% spectrum of Ass
[S,d] = projectors(Ass);
% create ST and make an array of eigenvalues
GmT = zeros(size(S));
Lm  = ones(n,n,length(d));
for IND = 1:length(d)
    GmT(:,:,IND) = S(:,:,IND)';
    Lm(:,:,IND) = d(IND)*ones(n,n);    
end
D = zeros(m,m,length(d));
for j = 1:length(d)
    Gk = S(:,:,j);
    Lk = d(j);
    if Ts ~= 0
        % discrete-time
        if Dss == 0
            FACTOR = ones(n,n,length(d))./(ones(n,n,length(d)) - Lk*Lm);
        else
            FACTOR = ((ones(n,n,length(d)) - Lk*ones(n,n,length(d)))+(ones(n,n,length(d)) - Lm))./...
                     ((ones(n,n,length(d)) - Lk*Lm).*(ones(n,n,length(d)) - Lk*ones(n,n,length(d))).*(ones(n,n,length(d)) - Lm)); % from Eq.30
        end
    else        
        % continuous-time
        if Dss == 0
            FACTOR = - ones(n,n,length(d))./(Lk*ones(n,n,length(d)) + Lm);
        else
            FACTOR = ((Lk*Lm)+Lk*ones(n,n,length(d)) + Lm)./((Lk*Lm).*(Lk*ones(n,n,length(d)) + Lm));
        end
    end
    D(:,:,j) = Css*(Gk*SIGMA)*sum(GmT.*FACTOR,3)*Css';
end
% zero lag covariance matrix
G = real(sum(D,3));
if abs(eigs(G,1)) < 0
    warning('mdac:dispersionAnalysis','Zero - lag covariance matrix is negative definite.')
end
% calculate percentage dispersions and other stuff
[wn,zn,Delta,normDeltaL2,normDeltaLinf,eigenVal] = getDA(d,D,G,Ts);
if nargout==0
    printDAonworkspace(wn,zn,normDeltaL2,normDeltaLinf,eigenVal,Ts);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------------------------------------------------------
function [wn,zn,Delta,normDeltaL2,normDeltaLinf,eigenVal] = getDA(d,D,G,Ts)

%
% returns,
%
%   1. natural frequencies (rad/s)
%   2. damping ratios
%   3. percentage modal dispersion matrices...
%   4. ... and their norms
%

% i. duplicate inputs and define data
wn = zeros(1,1);  zn = wn;  eigenVal = wn;
Delta = zeros(size(D,1),size(D,2),1); normDeltaL2 = wn; normDeltaLinf = wn;
k = 1;
ind = 1;
% ii. calculate
while k <= length(d)
    if isreal(d(k)) == 1
        wn(ind,1) = 0;
        zn(ind,1) = 1;
        Delta(:,:,ind) = 100*D(:,:,k)./G;
        eigenVal(ind,1) = d(k);
        k = k+1;
    else        
        eigenVal(ind,1) = d(k);
        if Ts ~= 0
            ctEig = log(d(k))/Ts;
        else
            ctEig = d(k);
        end
        wn(ind,1) = sqrt(real(ctEig)^2 + imag(ctEig)^2);
        zn(ind,1) = -real(ctEig)/wn(ind,1);        
        Delta(:,:,ind) = 100*(real(D(:,:,k)+D(:,:,k+1)))./G;
        k = k+2;
    end
    normDeltaL2(ind,1) = norm(Delta(:,:,ind),2);%max(max(abs(Delta(:,:,ind))));
    normDeltaLinf(ind,1) = norm(Delta(:,:,ind),inf);
    ind = ind+1;
end
% sort
[MAT,ind] = sortrows([wn/(2*pi) zn normDeltaL2 normDeltaLinf eigenVal],1);
wn = MAT(:,1);
zn = MAT(:,2);
normDeltaL2 = MAT(:,3);
normDeltaLinf = MAT(:,4);
eigenVal = MAT(:,5);
Delta = Delta(:,:,ind);

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
function printDAonworkspace(wn,zn,normDeltaL2,normDeltaLinf,eigenVal,Ts)

%
% for no output arguments, the
% function prints the results
% on workspace.
%

if Ts == 0
    STR = 'Eigenvalue';
else
    STR = 'Discrete-time eigenvalue';
end
form = '%4.3f';
% dispersion rate
deltaRateL2 = normDeltaL2/max(normDeltaL2);
deltaRateLinf = normDeltaLinf/max(normDeltaLinf);
% sort
MAT = [wn zn deltaRateL2 deltaRateLinf eigenVal];
% Print on screen via DPRINT    
mNo = dprint(transpose(linspace(1,length(wn),length(wn))),'Mode No.','%02g');
wP = dprint(MAT(:,1),'Frequency (Hz)',form);
zP = dprint(100*MAT(:,2),'Damping (%)',form);
DeltaP2 = dprint(MAT(:,3),'L2 Dispersion',form);
DeltaPinf = dprint(MAT(:,4),'Linf Dispersion',form);
poleReal = num2str(real(MAT(:,5)),form);
poleImag = num2str(abs(imag(MAT(:,5))),form);
poleStruc = blanks(size(poleReal,2)+size(poleImag,2)+6);
for k = 1:size(poleReal,1)
    % create structured vector of eigenvalues
    poleStruc(k,:) = [poleReal(k,:),' +/- ',poleImag(k,:),'i'];
end
eigP = dprint(poleStruc,STR,'%4.3f');
disp([mNo wP zP DeltaP2 DeltaPinf eigP])
%--------------------------------------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




























    
