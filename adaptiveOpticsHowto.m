%% ADAPTIVE OPTICS M0DELING WITH OOMAO
% Demonstrate how to build a simple adaptive optics system

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
    'fieldOfViewInArcMin',2.5,...
    'resolution',nPx,...
    'samplingTime',1/100);

%% Wavefront sensor
nLenslet = 10;
wfs = shackHartmann(nLenslet,nPx,0.75);
setValidLenslet(wfs,utilities.piston(nPx))

%% Source
ngs = source;

%% Deformable mirror
nActuator = nLenslet + 1;
bif = influenceFunction('monotonic',25/100);
dm = deformableMirror(nActuator,...
    'modes',bif,...
    'resolution',nPx,...
    'validActuator',wfs.validActuator);


%% Building the system
ngs=ngs.*tel*dm*wfs;
wfs.referenceSlopes = wfs.slopes;
grabAndProcess(wfs)
slopesAndFrameDisplay(wfs)
wfs.slopesListener.Enabled = true;
wfs.camera.frameListener.Enabled = true;

% tel.opticalAberration = atm;

%% Browsing the DM
wfs.slopesListener.Enabled = true;
wfs.camera.frameListener.Enabled = true;
wfs.paceMaker.Period = 0.1;
% start(wfs.paceMaker)

figure
% slopesDisplay(wfs,'parent',subplot(1,2,1))
imagesc(dm)%,'parent',subplot(1,2,2))
dm.surfaceListener.Enabled = true;
for k=1:97;
    dm.coefs(k)=10;
    ngs=ngs.*tel*dm*wfs;
%     pause(0.5);
    drawnow
    dm.coefs(k)=0;
end

% stop(wfs.paceMaker)

%% DM/WFS calibration
wfs.slopesListener.Enabled = false;
wfs.camera.frameListener.Enabled = false;
dm.surfaceListener.Enabled = false;
dm.coefsDefault = 0;
stroke = 3;
dm.coefs = eye(dm.nValidActuator)*stroke;
ngs=ngs.*tel*dm*wfs;
calibrationMatrix = wfs.slopes./stroke;
figure(10)
subplot(1,2,1)
imagesc(calibrationMatrix)
xlabel('DM actuators')
ylabel('WFS slopes [px]')
ylabel(colorbar,'slopes/actuator stroke')

% %% OLD library
% wfs = shackHartmann(nLenslet,0.8,6);
% cif = clampedInfluenceFunction(25/100,2/nLenslet);
% [x,y,r,o] = cartAndPol(nPx);
% telPupil = piston(nPx);
% wfs = data(wfs,telPupil,'reference');
% dm = zonalDeformableMirror(nActuator,cif,x,y);
% dm = set(dm,'pupilFootprint',telPupil,'validActuator',get(wfs,'validActuator'));
% dm = set(dm,'validActuator',get(wfs,'validActuator'));
% calib = calibration(wfs,dm);
% figure
% imagesc(calib)
% 
% wfs = data(wfs,telPupil);
% r = r./max(r(logical(telPupil)));
% zern = zernikePolynomials(1:6,r,o);
% zern = set(zern,'coefficients',ones(6,1));
% wfs = data(wfs,get(zern,'wave2D'));
% z = zernikeCoefsEstimate(wfs);

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

%% closed loop
wfs.slopesListener.Enabled = false;
wfs.camera.frameListener.Enabled = false;
gain = 0.5;
tel = tel+atm;
dm.coefs = zeros(dm.nValidActuator,1);
ngs=ngs.*tel;
turbPhase = ngs.meanRmPhase;
ngs=ngs*dm*wfs;
figure(11)
h = imagesc([turbPhase,ngs.meanRmPhase,dm.phase]);
axis equal tight
colorbar
pause
ngs.bufSeq = true(1,2);
while true
    ngs=ngs.*+tel; 
    turbPhase = ngs.meanRmPhase;
    ngs=ngs*dm*wfs; 
    residualDmCoefs = commandMatrix*wfs.slopes;
    dm.coefs = dm.coefs - gain*residualDmCoefs;
    set(h,'Cdata',[turbPhase,ngs.meanRmPhase,-dm.phase])
    
    drawnow
end

%% Zernike measurement
maxRadialDegree = 8;
zern = zernike(2:zernike.nModeFromRadialOrder(maxRadialDegree),'resolution',nPx,'pupil',tel.pupil);
zern.lex = false;
% figure(10)
% imagesc(zern.phase)
zern.c = eye(zern.nMode);
% wfs.lenslets.wave = zern.wave;
% grabAndProcess(wfs)
ngs=ngs.*zern*wfs;
% slopesAndFrameDisplay(wfs)
% z = getZernike(wfs,maxRadialDegree);
z = zernike(1:zernike.nModeFromRadialOrder(maxRadialDegree))\wfs;
Dz = z.c(2:end,:);

%% With noise
onAxis = source;
onAxis.wavelength = photometry.R;
onAxis.magnitude = 16;
% wfs.lenslets.lightSource = onAxis;
wfs.camera.readOutNoise = 10;
wfs.camera.photonNoiseLess = false;
tel=tel-atm;
% wfs.lenslets.wave = tel;
wfs.framePixelThreshold = 0;
% grabAndProcess(wfs)
onAxis=onAxis.*tel*wfs;
slopesAndFrameDisplay(wfs)

%% noise convariance matrix
nMeas = 1000;
slopes = zeros(wfs.nSlope,nMeas);
for kMeas=1:nMeas
%     grabAndProcess(wfs)
    onAxis=onAxis.*tel*wfs;
    slopes(:,kMeas) = wfs.slopes;
end
Cn = slopes*slopes'/nMeas;
figure(5)
subplot(1,2,1)
imagesc(Cn)
axis equal tight
colorbar
wfs.slopes = slopes;
% z = getZernike(wfs,maxRadialDegree);
z = z\wfs;
Czn = z.c(2:end,:)*z.c(2:end,:)'/nMeas;
subplot(1,2,2)
imagesc(Czn)
axis equal tight
colorbar

% %% Phase reconstruction (BROKEN)
% tel = tel+atm;
% %% wavefront reconstruction least square fit
% offAxis = reset(onAxis);%source('zenith',0*cougarConstants.arcmin2radian,'azimuth',0);
% offAxis=offAxis.*tel;
% ps = offAxis.meanRmPhase;
% % wfs.lenslets.wave     = tel;
% % wfs.camera.readOutNoise = 1;
% % grabAndProcess(wfs)
% offAxis=offAxis*wfs;
% % z = getZernike(wfs,maxRadialDegree);
% z = z\wfs;
% zern.c = Dz\z.c(2:end);
% ngs = source.*zern;
% phaseLS = ngs.phase;
% 
% %% wavefront reconstruction minimum variance
% Cz = phaseStats.zernikeCovariance(zern,atm);
% % M = Cz*Dz'/(Dz*Cz*Dz'+Czn);
% M = Cz*Dz'/(Dz*Cz*Dz'+Czn);
% zern.c = M*z.c(2:end);
% ngs = source.*zern;
% phaseMV = ngs.phase;
% figure(11)
% subplot(2,1,1)
% imagesc([ps,phaseLS,phaseMV])
% axis equal tight
% colorbar
% subplot(2,1,2)
% imagesc([ps-phaseLS,ps-phaseMV])
% axis equal tight
% colorbar
% 
% 
