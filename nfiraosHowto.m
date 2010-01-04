%% NFIRAOS HOWTO
% How to define object to simulate NFIRAOS

%% The telescope
tmt = telescope(30,...
    'fieldOfViewInArcMin',2,...
    'resolution',60*15,...
    'samplingTime',1/800);
% tmt = telescope(30,'resolution',60*15);
tmt.focalDistance = 90e3;

%% A lenslet array
la = lensletArray(60);
la.nyquistSampling = 0.5;
la.wave = tmt;
la.imageletsListener.Enabled = true;
propagateThrough(la)
% imagesc(la)

%% One LGS on axis
lgs = source('height',linspace(86,94,15).*1e3,'wavelength',589e-9);
la.lightSource = lgs;
propagateThrough(la)

%% 6 LGSs
lgs = source('asterism',{[0,0],[5,60*cougarConstants.arcsec2radian,0]},...
    'height',linspace(86,94,15).*1e3,'wavelength',589e-9);
la.lightSource = lgs;
propagateThrough(la)

%% LGS + phase screens
atm = atmosphere(0.65e-6,0.15,60,...
    'altitude',2e3,...
    'fractionnalR0',1,...
    'windSpeed',10,...
    'windDirection',0);
tmt.opticalAberration = atm;
propagateThrough(la)

