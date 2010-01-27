%% ADAPTIVE OPTICS TOMOGRAPHY HOWTO
% Demonstrate how to build a tomographic adaptive optics system

%% Atmosphere 
Cn2 = [6.39 3.94 1.46 1.73 3.11 2.69 2.81];
fr0 = Cn2/sum(Cn2);
atm = atmosphere(photometry.V,0.2,30,...
    'altitude',[0,0.5,1,2,4,8,16]*1e3,...
    'fractionnalR0',fr0,...
    'windSpeed',[5.6,5.8,6.3,7.6,13.3,19.1,12.1],...
    'windDirection',[0,5,15,30,60,90,180]*pi/180);
atm.wavelength = photometry.R;

%% Telescope
nPx = 8*16;
tel = telescope(8.2,...
    'fieldOfViewInArcMin',2.5,...
    'resolution',nPx,...
    'samplingTime',1/500);


%% Wavefront sensor
nLenslet = 16;
wfs = shackHartmann(nLenslet,nPx,0.75);
setValidLenslet(wfs,utilities.piston(nPx))

%% Deformable mirror
nActuator = nLenslet + 1;
bif = influenceFunction('monotonic',25/100);
dm = deformableMirror(nActuator,...
    'modes',bif,...
    'resolution',nPx,...
    'validActuator',wfs.validActuator);

%% Zernike definition 
maxRadialDegree = 12;
zernModeMax = zernike.nModeFromRadialOrder(maxRadialDegree);
zern = zernike(2:zernModeMax,tel.D,'resolution',nPx);

%% Sources
 % Calibration source
ngs = source('wavelength',photometry.R);
gs = source('asterism',{[3,0.5*constants.arcmin2radian,0]},...
    'wavelength',photometry.R,'magnitude',0);
% scs = source('asterism',{[0,0],[30*cougarConstants.arcsec2radian,pi/4]},'wavelength',photometry.H);
scs = source('asterism',{[0,0],...
    [12, 30*cougarConstants.arcsec2radian, 0],...
    [12, 60*cougarConstants.arcsec2radian, 0]},'wavelength',photometry.H);
nScs = length(scs);
nGs = length(gs);

%% Building the system
ngs=ngs.*tel*dm*wfs;
wfs.referenceSlopes = wfs.slopes;
slopesAndFrameDisplay(wfs)
% pause
+ngs;
slopesAndFrameDisplay(wfs)

%   %% DM/WFS calibration
% dm.coefsDefault = 0;
% stroke = 3;
% dm.coefs = eye(dm.nValidActuator)*stroke;
% +ngs;
% dm.coefs = 0;
% calibrationMatrix = wfs.slopes./stroke;
% figure(10)
% subplot(1,2,1)
% imagesc(calibrationMatrix)
% xlabel('DM actuators')
% ylabel('WFS slopes [px]')
% ylabel(colorbar,'slopes/actuator stroke')
% 
% %% Command matrix derivation
% [nS,nC] = size(calibrationMatrix);
% [U,S,V] = svd(calibrationMatrix);
% nThresholded = 4;
% eigenValues = diag(S);
% subplot(1,2,2)
% semilogy(eigenValues,'.')
% xlabel('Eigen modes')
% ylabel('Eigen values')
% iS = diag(1./eigenValues(1:end-nThresholded));
% iS(nC,nS) = 0;
% commandMatrix = V*iS*U';

%% Zernike measurement
zern.lex = false;
% figure(10)
% imagesc(zern.phase)
zern.c = eye(zern.nMode);
ngs=ngs.*zern*wfs;
z = zernike(1:zernike.nModeFromRadialOrder(maxRadialDegree),tel.D)\wfs;
Dz = z.c;

%% With noise
ngs.magnitude = 0;
wfs.camera.readOutNoise = 0;
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
Czn = z.c*z.c'/nMeas;
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
figure(3)
subplot(1,atm.nLayer+1,[1,atm.nLayer])
imagesc(tel,[gs,scs])
subplot(1,atm.nLayer+1,1+atm.nLayer)
polar(scs)
hold on
polar(gs,'ro')
hold off
% Zernike section expansion
%% Turbulence covariance matrix
tic
S = phaseStats.zernikeAngularCovariance(zern,atm,gs);
toc
S = cell2mat(S);
%% Data/Target covariance
C = phaseStats.zernikeAngularCovariance(zern,atm,gs,scs);
%% tomographic matrices
CznAst = blkdiag( Czn , Czn , Czn );
DzAst = blkdiag( Dz , Dz , Dz );

% %% wavefront reconstruction least square fit
% gs=gs.*tel;
% ps = [gs.meanRmPhase];
% gs=gs*wfs
% z = z\wfs;
% zern.c = Dz\z.c(2:end,:);
% phaseLS = reshape(zern.phase,nPx,nPx*length(gs));
% 
% %% wavefront reconstruction minimum variance
% M = S/(S+CznAst);
% M = S*DzAst'/(DzAst*S*DzAst'+CznAst);
% z = z - 1; % piston removed
% zern.c = reshape(M*z.c(:),z.nMode,[]);
% phaseMV = reshape(zern.phase,nPx,nPx*length(gs));
% figure(12)
% imagesc([ps;phaseLS;phaseMV])
% axis equal tight xy

%% Target matrix
M = cell(nScs,1);
denom = DzAst*S*DzAst'+CznAst;
for kScs = 1:nScs
    M{kScs,1} = cell2mat(C(:,kScs))'*DzAst'/denom;
end
M = cell2mat(M);

%% Command matrix
lambdaRatio = gs(1).wavelength/scs(1).wavelength;
gs=gs.*tel*wfs;
z = zernike(1:zernModeMax)\wfs;
zern.c = reshape(M*(lambdaRatio*z.c(:)),z.nMode,nScs);
ngs = ngs.*zern;
scs = scs.*tel;
turbPhase = [scs.meanRmPhase];
nIt =100;
turbPhaseStd = zeros(nIt,nScs);
turbPhaseStd(1,:) = scs.var;
figure
imagesc([scs.meanRmPhase;reshape(ngs.phase,nPx,[])])
axis equal tight xy
colorbar

%% DM
zern.lex = true;
zern2dm = dm.modes.modes(tel.pupilLogical,:)\zern.p(tel.pupilLogical,:)/2;%\zPoly/2;
dm.coefs = zern2dm*zern.c;

turbRes = zeros(nPx,nPx*nScs,nIt);
turbRes(:,:,1) = turbPhase;
turbResStd = turbPhaseStd;
figure
% plot([turbPhaseStd(1:k,:),turbResStd(1:k,:)],'.');
% set(h(1),'YDataSource',turbPhaseStd(:,1))
h = imagesc([turbPhase;turbRes(:,:,1);lambdaRatio*reshape(-dm.phase,nPx,[])]);
axis equal tight xy
colorbar

%% Open-Loop
k = 1;
scs(1).saveImage = true;
scs(2).saveImage = true;
log = logBook.checkIn;
log.verbose=false;
warning off MATLAB:rankDeficientMatrix
tic
while k<nIt
    
    % propagation of the guide stars to the WFSs
    gs=gs.*tel*wfs;
    % slopes projection onto the Zernikes, piston removed
%     z = z\wfs;
%     z.c = wfs.zernCoefs;
    % DMs command vectors
    dm.coefs = lambdaRatio*zern2dm*reshape(M*wfs.zernCoefs(:),z.nMode,nScs);
    % Atmosphere update
    +tel;
    % propagation of science star to the telescope through the atmosphere
    scs = scs.*tel;
    turbPhase = [scs.meanRmPhase];
    k = k + 1 ;
    turbPhaseStd(k,:) = scs.var;
    % propagation of science star resumes to the DMs
    scs = scs*dm;
    turbRes(:,:,k) = [scs.meanRmPhase];
    
    turbResStd(k,:) = scs.var;
%     set(h,'Ydata',[turbPhaseStd,turbResStd])
%     set(h,'Cdata',[turbPhase;turbRes(:,:,k);reshape(-dm.phase,nPx,[])])
%     drawnow
    
end
toc

%%
u = [1 nIt];
atmWavelength = atm.wavelength;
atm.wavelength = scs(1).wavelength;
figure
plot(turbPhaseStd,'.')
hold on
plot(turbResStd,'.')
hold off
line(u,phaseStats.zernikeResidualVariance(1,atm,tel)*ones(2,1),'color','k','LineWidth',2)
line(u,ones(2,1)*mean(turbPhaseStd),'lineStyle','--')
line(u,phaseStats.zernikeResidualVariance(zern.nMode,atm,tel)*ones(2,1),'color','r')
grid
xlabel('Iterations')
ylabel('Variance [rd^2]')
atm.wavelength = atmWavelength;

%%
turbRes = reshape(turbRes,nPx,nPx*nScs*nIt);
turbRes = mat2cell(turbRes,nPx,nPx*ones(1,nScs*nIt));
%%
residualWave = cellfun(@(x)tel.pupil.*exp(1i*x),turbRes,'UniformOutput',false);
% Optical transfer function
normTel = sum(tel.pupil(:));
nOtf = 2*nPx;
otfTel = fftshift(ifft2(abs(fft2(tel.pupil,nOtf,nOtf)).^2))/normTel;
bigRam = false;
if bigRam
    tic
    otfPd = cellfun(@(x)fftshift(ifft2(abs(fft2(x,nOtf,nOtf)).^2))/normgTel,residualWave,'UniformOutput',false);
    toc
    clear residualWave
    for kScs = 1:nScs
        otfPd{kScs} = mean(reshape(cell2mat(otfPd((nScs+kScs):nScs:end)),nOtf,nOtf,nIt-1),3);
    end
    meanOtfPd = otfPd(1:nScs);
    clear otfPd
else
    meanOtfPd = cell(1,nScs);
    h = waitbar(0,'Computing OTFs ...');
    tic
    for kScs = 1:nScs
        otfPd = cellfun(@(x)fftshift(ifft2(abs(fft2(x,nOtf,nOtf)).^2))/normTel,...
            residualWave((nScs+kScs):nScs:end),'UniformOutput',false);
        meanOtfPd{kScs} = mean(reshape(cell2mat(otfPd),nOtf,nOtf,nIt-1),3);
        waitbar(kScs/nScs)
    end
    toc
    close(h)
    clear residualWave otfPd
end
% Strehl ratio
u = linspace(-tel.D,tel.D,nOtf);
strehlRatioFun = @(x)real(trapz(u,trapz(u,x)))/tel.area;
strehlRatio = cellfun(strehlRatioFun,meanOtfPd);
% entraped energy
[x,y] = meshgrid(u);
ircsSlitWidthInArcsec = 0.14';
a = (ircsSlitWidthInArcsec/(photometry.H/tel.D*constants.radian2arcsec))/tel.D; % diameter
eeFilter = a^2*(sin(pi.*x.*a)./(pi.*x.*a)).*...
    (sin(pi.*y.*a)./(pi.*y.*a));
eNrgFun = @(x) real(trapz(u,trapz(u,x.*eeFilter)));
fprintf(' > Diffraction limited nrg (fine/coarse) : %4.1f%%/%4.1f%%\n',...
    entrappedEnergy(tel-atm,a/2,'square','otf')*100,eNrgFun(otfTel)*100)
tel = tel + atm;
eNrg = cellfun(eNrgFun,meanOtfPd);

%%
[x,y] = pol2cart([scs.azimuth],[scs.zenith]*cougarConstants.radian2arcsec);
z = zeros(size(x));
tri = delaunay(x,y);
figure
subplot(1,2,1)
trisurf(tri,x,y,z,strehlRatio)
view(2)
shading interp
axis square
colorbar
hold on
polar(scs,'k*')
polar(gs,'wo')
hold off
title('Strehl ratio')
set(gca,'View',[0 90],'Box','on')
subplot(1,2,2)
trisurf(tri,x,y,z,eNrg)
view(2)
shading interp
axis square
colorbar
hold on
polar(scs,'k*')
polar(gs,'wo')
hold off
title('entrapped energy')
set(gca,'View',[0 90],'Box','on')
% [o,r] = meshgrid([0:12]*pi/6,[0 30 60]);
% [x,y] = pol2cart(o,r);
% u = reshape(strehlRatio(2:end),[],2)';
% u = [repmat(strehlRatio(1),1,13); u, u(:,1)];
% v = reshape(eNrg(2:end),[],2)';
% v = [repmat(eNrg(1),1,13); v, v(:,1)];
% figure
% subplot(1,2,1)
% pcolor(x,y,u)
% shading interp
% axis square
% colorbar
% subplot(1,2,2)
% pcolor(x,y,v)
% shading interp
% axis square
% colorbar
