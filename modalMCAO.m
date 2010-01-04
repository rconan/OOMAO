%% Modal MCAO
% % atmosphere definition: wavelength,r0,L0,...
% atm   = atmosphere(2.2e-6,0.8,25,...
%     'altitude',[0.       1800. 3300. 5800. 7400. 13100. 15800],...
%     'fractionnalR0',[0.646   0.078  0.119 0.035 0.025  0.080 0.017]);
% nRadialOrder = 2;
% % dms definition:
% % {firstMode:lastMode},{telescopeDiameter,fieldOfView},conjugationAltitude
% dm(1) = zernikeDeformableMirror({2:zernike.nModeFromRadialOrder(nRadialOrder)},{8,'fieldOfViewInArcsec',30});
% dm(2) = zernikeDeformableMirror({4:zernike.nModeFromRadialOrder(nRadialOrder)},{8,'fieldOfViewInArcsec',30},10e3);
% % wfs guide start: at Matlab prompt: help source
% gs    = source('asterism',{[3,15.*cougarConstants.arcsec2radian,0]});
% % wfs definition: firstMode:lastMode 
% wfs   = zernike(2:zernike.nModeFromRadialOrder(nRadialOrder));
% % optimization direction: at Matlab prompt: help source
% pd    = source('asterism',{[0,0],[3,15.*cougarConstants.arcsec2radian,0]});

% % NFIRAOS
% atmosphere definition: wavelength,r0,L0,...
fovInArcsec = 30;
atm   = atmosphere(2.2e-6,0.8,25,...
    'altitude',[0.       1800. 3300. 5800. 7400. 13100. 15800],...
    'fractionnalR0',[0.646   0.078  0.119 0.035 0.025  0.080 0.017]);
nRadialOrder = 2;
% dms definition:
% {firstMode:lastMode},{telescopeDiameter,fieldOfView},conjugationAltitude
dm(1) = zernikeDeformableMirror({2:zernike.nModeFromRadialOrder(nRadialOrder)},{8,'fieldOfViewInArcsec',fovInArcsec});
dm(2) = zernikeDeformableMirror({4:zernike.nModeFromRadialOrder(nRadialOrder)},{8,'fieldOfViewInArcsec',fovInArcsec},10e3);
% wfs guide start: at Matlab prompt: help source
gs    = source('asterism',{[5,0.5*fovInArcsec.*cougarConstants.arcsec2radian,0]});
% wfs definition: firstMode:lastMode 
wfs   = zernike(2:zernike.nModeFromRadialOrder(nRadialOrder));
% optimization direction: at Matlab prompt: help source
pd    = source('asterism',{[0,0],[5,0.5*fovInArcsec*cougarConstants.arcsec2radian,0]});


nDm  = length(dm);
nGs  = length(gs);
nPd  = length(pd);
% optimization direction weight
wpd  = ones(1,nPd)/nPd;

fprintf('-------- MODAL MCAO STARTED --------\n')
tic

%% Projection
P = cell(nPd,nDm);
o = linspace(0,2*pi,101);
x0 = cos(o);
y0 = sin(o);
figure
for kDm = 1:nDm
    fprintf('@(Projection)> ')
    subplot(1,nDm,kDm)
    title(sprintf('DM #%d',kDm))
    line(x0,y0)
    for kPd = 1:nPd
        %         D = tel.diameterAt(dm(kDm).height);
        fprintf('pd#%d/dm#%d - ',kPd,kDm)
        P{kPd,kDm} = footprintProjection(dm(kDm),pd(kPd));
        [P{kPd,kDm},x,y] = footprintProjection(dm(kDm),pd(kPd));
        line(x,y,'color','r','linestyle','--')
    end
    axis square
    fprintf('\b\b\b\n')
end

%% Data covariance
S = cell(nGs);
for iGs = 1:nGs
    fprintf('@(Data covariance)> ');
    gsCurrent = gs(iGs);
    parfor jGs = iGs:nGs
        fprintf('gs#%d/gs#%d - ',iGs,jGs)
        S{iGs,jGs} = phaseStats.zernikeAngularCovariance(wfs,atm,dm,[gsCurrent,gs(jGs)]);
    end
    fprintf('\b\b\b\n')
end
index = cellfun(@isempty,S);
S(index) = cellfun(@transpose,S(triu(~index,1)),'UniformOutput',false);
S = cell2mat(S);

%% Data/Target covariance
C = cell(nGs,nPd);
for kPd=1:nPd
    fprintf('@(Data/Target covariance)> ');
    pdCurrent= pd(kPd);
    parfor kGs=1:nGs
        fprintf('pd#%d/gs#%d - ',kPd,kGs)
        C{kGs,kPd} = phaseStats.zernikeAngularCovariance(wfs,atm,dm,[gs(kGs),pdCurrent]);
    end
    fprintf('\b\b\b\n')
end

%% Target matrix
T = 0;
R = 0;
for kPd = 1:nPd
    PnDM = cell2mat(P(kPd,:));
    T = T + wpd(kPd)*cell2mat(C(:,kPd))*PnDM(2:end,:);
    R = R + wpd(kPd)*(PnDM'*PnDM);
end

toc
fprintf('--------- MODAL MCAO ENDED ---------\n')

%% Command matrix
% M = pinv(R,1e-6)*T'/S;
M = R\T'/S;

%% Results
pistonFreePhaseVariance = phaseStats.variance(atm) - phaseStats.zernikeVariance(zernike(1),atm,dm);
scaoVariance = phaseStats.variance(atm) - ...
    sum(phaseStats.zernikeVariance(zernike(1:zernike.nModeFromRadialOrder(nRadialOrder)),atm,dm));
mcaoVariance = pistonFreePhaseVariance - trace(2*M*T - R*M*S*M');
mcaoTargetVariance = zeros(1,nPd);
for kPd = 1:nPd
    PnDM = cell2mat(P(kPd,1:nDm));
    T = cell2mat(C(1:nGs,kPd))*PnDM(2:end,:);
    R = (PnDM'*PnDM);
    mcaoTargetVariance(kPd) = pistonFreePhaseVariance - trace(2*M*T - R*M*S*M');    
end
fprintf('--------- MODAL MCAO RESULTS ---------\n')
fprintf('SCAO variance: %3.2f rd^2\nMCAO variance: %3.2f rd^2\n',scaoVariance,mcaoVariance)
for kPd = 1:nPd
    fprintf('MCAO variance target#%d: %3.2f rd^2\n',kPd,mcaoTargetVariance(kPd))
end
