%% ADAPTIVE OPTICS MODELING WITH OOMAO - adaptation to test the integration of the pyramid wave fromt sensor
% Demonstrate how to build a simple closed-loop single conjugated adaptive
% optics system

%%
%close all
clear all
clc 
close all

%%
atm = atmosphere(photometry.V,.15,30,...
    'altitude',0,...
    'fractionnalR0',1,...
    'windSpeed',100,...
    'windDirection',0);

%% Definition of the telescope
nPx = 40;
tel = telescope(3.6,...
    'fieldOfViewInArcMin',2.5,...
    'resolution',nPx,...
    'samplingTime',1/100);

%% Definition of a calibration source
ngs = source('wavelength',photometry.J); 

ngs_sh = source('wavelength',photometry.J);
%ngs_sh = source('asterism',{[2,arcsec(45),0]},'magnitude',10); 

%% SH
nLenslet = 10;
wfs_sh = shackHartmann(nLenslet,nPx,0.75);
wfs_sh.camera.readOutNoise = 1;


%% Definition of the wavefront sensor
%Experimental Pyramid WFS, expect some rough edges and maybe some bugs
%The pyramid takes only one argument, which is the pixel resolution of the
%telescope it is associated with.


pyr = pyramid(nLenslet,nPx,'modulation',2,'binning',1);
pyr.camera.readOutNoise = 0;



%%
% Propagation of the calibration source to the WFS through the telescope
ngs = ngs.*tel*pyr;
ngs_sh = ngs_sh.*tel*wfs_sh;

%%
%Calibration of the sensor on the current light
%wfs.setmodulation(3)
pyr.INIT
wfs_sh.INIT

%%
% A new frame read-out and slopes computing:
+pyr;
+wfs_sh;
%%
% The WFS camera display:
figure(1)
imagesc(pyr.camera)
figure(2)
imagesc(wfs_sh.camera)
%%
% The WFS slopes display:
figure(3)
slopesDisplay(pyr)
figure(4)
slopesDisplay(wfs_sh)


%%
%Propagating through the atmosphere
% tel = tel+atm;
% 
% for ii = 1 :10
% ngs=ngs.*+tel;
% end
% ngs_sh=ngs_sh.*tel;
% 
% ngs=ngs*wfs;
% ngs_sh=ngs_sh*wfs_sh;
% 
% % A new frame read-out and slopes computing:
% +wfs;
% +wfs_sh;

%% 
zer = zernike(3,tel.D, 'resolution', nPx);
zer.c = 0.1/ngs.waveNumber;
ngs = ngs.*tel*zer*pyr;
ngs_sh = ngs_sh.*tel*zer*wfs_sh;

% A new frame read-out and slopes computing:
+pyr;
+wfs_sh;

%%
% The WFS camera display:
figure(1)
imagesc(pyr.camera)
figure(2)
imagesc(wfs_sh.camera)
%%
% The WFS slopes display:
figure(3)
slopesDisplay(pyr)
figure(4)
slopesDisplay(wfs_sh)
% 
% wfs.hilbertSlopes(tel,ngs)
% figure(5)
% slopesDisplay(wfs)

%% Gain Calibration
for i = 1:20
    zer.c = (i-10)*0.5/ngs.waveNumber;
    ngs = ngs.*tel*zer*pyr;
    sx(i) = mean(pyr.slopes(1:end/2));
    sy(i) = mean(pyr.slopes(end/2+1:end));
end
syTh = 4*([1:20]-10)*0.5;%/ngs.waveNumber;
figure,hold
plot(syTh, sy,'o')
plot(syTh, syTh)

%% Poke matrix
bifa = influenceFunction('monotonic',0.4);

dm = deformableMirror(nLenslet+1,'modes',bifa,...
    'resolution',tel.resolution,...
    'validActuator',wfs_sh.validActuator);


%ngs = ngs.*tel;
%calibDm = calibration(dm,pyr,ngs,ngs.wavelength,nLenslet+1,'cond',1e2);


pyr.binning = 2;
ngs = ngs.*tel*pyr;
pyr.INIT


% brute-force approach
ngs = ngs.*tel*dm*pyr;
for iAct = 1:dm.nValidActuator
    c = zeros(dm.nValidActuator,1);
    c(iAct) = 1e-6;
    dm.coefs = c;
    +ngs;
    pokeMatrix(iAct,:) = pyr.slopes;
end

% OOMAO-compatible 
dm.coefs = eye(dm.nValidActuator)*ngs.wavelength/10;
ngs.*tel*dm;
calib = calibration(dm,pyr,ngs,ngs.wavelength/10,90);



% modal interaction matrix
zer = zernike(2:500,tel.D, 'resolution', nPx);
ngs = ngs.*tel*zer*pyr;
pokeModalMatrix = [];
for kZer=1:zer.nMode
    kZer
    c = zeros(zer.nMode,1);
    c(kZer) = 10e-9;
    zer.c = c;
    tic
    +ngs;
    toc
    pokeModalMatrix(:,kZer) = pyr.slopes;
end
figure
loglog(diag(pokeModalMatrix'*pokeModalMatrix))
hold on

% SH modal interaction matrix
zer.c = eye(zer.nMode)*1e-10;%/ngs.waveNumber;
tic,
ngs = ngs.*tel*zer*wfs_sh;
toc
pokeModalMatrixSh = wfs_sh.slopes;
loglog(diag(pokeModalMatrixSh'*pokeModalMatrixSh),'r')
hold on

%ipokeModalMatrixSh = pinv(pokeModalMatrixSh);
%loglog(diag(ipokeModalMatrixSh*ipokeModalMatrixSh'),'k--')
