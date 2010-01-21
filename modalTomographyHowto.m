%% ADAPTIVE OPTICS TOMOGRAPHY HOWTO
% Demonstrate how to build a tomographic adaptive optics system

%% Atmosphere 
atm = atmosphere(photometry.V,0.15,30,...
    'altitude',[0,4,10]*1e3,...
    'fractionnalR0',[0.7,0.25,0.05],...
    'windSpeed',[5,10,20],...
    'windDirection',[0,pi/4,pi]);
atm.wavelength = photometry.R;

%% Telescope
nPx = 60;
tel = telescope(3.6,...
    'fieldOfViewInArcMin',1,...
    'resolution',nPx,...
    'samplingTime',1/100);

%% Sources
ngs = source;
ast = source('asterism',{[3,30*cougarConstants.arcsec2radian,0]});

%% Wavefront sensor
nLenslet = 10;
wfs = shackHartmann(nLenslet,nPx,0.75);
setValidLenslet(wfs,utilities.piston(nPx))

%% Deformable mirror
nActuator = nLenslet + 1;
bif = influenceFunction('monotonic',25/100);
dm = deformableMirror(nActuator,...
    'modes',bif,...
    'resolution',nPx,...
    'validActuator',wfs.validActuator,...
    'zLocation',4e3);

%% Building the system
ngs=ngs.*tel*dm*wfs;
wfs.referenceSlopes = wfs.slopes;
slopesAndFrameDisplay(wfs)
% pause
+ngs;
slopesAndFrameDisplay(wfs)

  %% DM/WFS calibration
dm.coefsDefault = 0;
stroke = 3;
dm.coefs = eye(dm.nValidActuator)*stroke;
+ngs;
dm.coefs = 0;
calibrationMatrix = wfs.slopes./stroke;
figure(10)
subplot(1,2,1)
imagesc(calibrationMatrix)
xlabel('DM actuators')
ylabel('WFS slopes [px]')
ylabel(colorbar,'slopes/actuator stroke')

%% Command matrix derivation
[nS,nC] = size(calibrationMatrix);
[U,S,V] = svd(calibrationMatrix);
nThresholded = 4;
eigenValues = diag(S);
subplot(1,2,2)
semilogy(eigenValues,'.')
xlabel('Eigen modes')
ylabel('Eigen values')
iS = diag(1./eigenValues(1:end-nThresholded));
iS(nC,nS) = 0;
commandMatrix = V*iS*U';

%% Zernike measurement
maxRadialDegree = 3;
zern = zernike(2:zernike.nModeFromRadialOrder(maxRadialDegree),'resolution',nPx,'pupil',tel.pupil);
zern.lex = false;
figure(10)
imagesc(zern.phase)
zern.c = eye(zern.nMode);
ngs=ngs.*zern*wfs;
z = zernike(1:zernike.nModeFromRadialOrder(maxRadialDegree))\wfs;
Dz = z.c(2:end,:);

%% With noise
ngs.wavelength = photometry.R;
ngs.magnitude = 0;
wfs.camera.readOutNoise = 5;
wfs.camera.photonNoiseLess = false;
wfs.framePixelThreshold = 0;
ngs=ngs.*tel*wfs;
slopesAndFrameDisplay(wfs)

%% noise convariance matrix
nMeas = 250;
slopes = zeros(wfs.nSlope,nMeas);
for kMeas=1:nMeas
    grabAndProcess(wfs)
    slopes(:,kMeas) = wfs.slopes;
end
Cn = slopes*slopes'/nMeas;
figure(5)
subplot(1,2,1)
imagesc(Cn)
axis equal tight
colorbar
wfs.slopes = slopes;
z = z\wfs;
Czn = z.c(2:end,:)*z.c(2:end,:)'/nMeas;
subplot(1,2,2)
imagesc(Czn)
axis equal tight
colorbar

%% Phase reconstruction
tel = tel+atm;
% %% wavefront reconstruction least square fit
% ngs=ngs.*tel;
% ps = ngs.meanRmPhase;
% ngs=ngs*wfs;
% z = z\wfs;
% zern.c = Dz\z.c(2:end);
% phaseLS = zern.phase;
% 
% %% wavefront reconstruction minimum variance
% Cz = phaseStats.zernikeCovariance(zern,atm,tel);
% M = Cz*Dz'/(Dz*Cz*Dz'+Czn);
% zern.c = M*z.c(2:end);
% phaseMV = zern.phase;
% figure(11)
% subplot(2,1,1)
% imagesc([ps,phaseLS,phaseMV])
% axis equal tight
% colorbar
% subplot(2,1,2)
% imagesc([ps-phaseLS,ps-phaseMV])
% axis equal tight
% colorbar

%% TOMOGRAPHY
 %% Tomography
pd = source('zenith',0*cougarConstants.arcsec2radian);
figure(3)
imagesc(tel,[ast,pd])
wpd = ones(size(pd));
gs = ast;
nGs = length(ast);
nPd = length(pd);
zernModeMax = zernike.nModeFromRadialOrder(maxRadialDegree);
wfsMaxRadialDegree = zernModeMax;
zernWfs = zernike(2:zernModeMax,'resolution',nPx,'pupil',tel.pupil);
tel = tel+atm;
% Zernike section expansion
%% Projection
nDm = length(dm);
P = cell(nPd,nDm);
for kDm = 1:nDm
    fprintf('@(Projection)> ')
    for kPd = 1:nPd
        fprintf('pd#%d/dm#%d - ',kPd,kDm)
        src = pd(kPd);
        delta = dm.zLocation.*tan(src.zenith).*...
            [cos(src.azimuth),sin(src.azimuth)];
        delta = delta*2/tel.diameterAt(dm.zLocation);
        alpha = tel.diameterAt(dm.zLocation)./tel.D;
        P{kPd,kDm} = smallFootprintExpansion(zernWfs,delta,alpha);
    end
    fprintf('\n')
end
% Turbulence covariance matrix
SFile = sprintf('S%d__18-Jan-2010.mat',maxRadialDegree);
if exist(SFile,'file')==2
    load(SFile)
else
    S = phaseStats.zernikeAngularCovariance(zernWfs,atm,tel,ast);
    SFile = sprintf(['S%d__',date,'.mat'],maxRadialDegree);
    save(SFile,'S','maxRadialDegree')
    fprintf(' S matrix saved to %s\n',SFile)
end

%% wavefront reconstruction least square fit
ast = source('asterism',{[3,30*cougarConstants.arcsec2radian,0]},...
    'wavelength',photometry.R,'magnitude',16)
ast=ast.*tel;
ps = [ast.meanRmPhase];
ast=ast*wfs
z = z\wfs;
zern.c = Dz\z.c(2:end,:);
phaseLS = reshape(zern.phase,nPx,nPx*length(ast));

%% wavefront reconstruction minimum variance
S = cell2mat(S);
CznAst = blkdiag( Czn , Czn , Czn );
DzAst = blkdiag( Dz , Dz , Dz );
M = S/(S+CznAst);
M = S*DzAst'/(DzAst*S*DzAst'+CznAst);
z = z - 1; % piston removed
zern.c = reshape(M*z.c(:),z.nMode,[]);
phaseMV = reshape(zern.phase,nPx,nPx*length(ast));
figure(12)
imagesc([ps;phaseLS;phaseMV])
axis equal tight xy

%% Data/Target covariance
C = cell(nGs,nPd);
for kPd=1:nPd
    fprintf('@(Data/Target covariance)> ');
    pdCurrent= pd(kPd);
    for kGs=1:nGs
        fprintf('pd#%d/gs#%d - ',kPd,kGs)
        C{kGs,kPd} = phaseStats.zernikeAngularCovariance(zernWfs,atm,tel,[gs(kGs),pdCurrent]);
    end
    fprintf('\b\b\b\n')
end
%% Target matrix
T = 0;
R = 0;
for kPd = 1:nPd
    PnDM = cell2mat(P(kPd,:));
    T = T + wpd(kPd)*DzAst*cell2mat(C(:,kPd))*PnDM(2:end,:);
    R = R + wpd(kPd)*(PnDM'*PnDM);
end

%% Command matrix
M = R\T'/(DzAst*S*DzAst'+CznAst);
ast=ast.*tel*wfs;
z = zernike(1:zernike.nModeFromRadialOrder(maxRadialDegree))\wfs;
z = z - 1;
zern.c = M*z.c(:);
pd = pd.*tel;
figure
imagesc([zern.phase,pd.meanRmPhase])
axis equal tight xy
colorbar
