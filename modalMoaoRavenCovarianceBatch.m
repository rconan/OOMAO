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
    [12, 60*cougarConstants.arcsec2radian, 0]});
nScs = length(scs);
nGs = length(gs);

% %% Turbulence covariance matrix
% S = phaseStats.zernikeAngularCovariance(zern,atm,gs);
% S = cell2mat(S);
%% Data/Target covariance
tic
C = phaseStats.zernikeAngularCovariance(zern,atm,gs,scs);
toc
save('C12','C')
