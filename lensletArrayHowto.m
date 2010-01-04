%% LENSLETARRAY HOWTO
% Demonstrate the use of the <matlab:doc('lensletArray') lensletArray> class

%% A simple 10x10 lenslet array
nlenslet = 10;
nPx = 60;
la = lensletArray(nlenslet)
% a circular pupil
la.wave = utilities.piston(nPx);
% propagated throught the lenslet array
propagateThrough(la)
% and displayed
imagesc(la)
% auto-update of the display
la.imageletsListener.Enabled = true;

%% changing the sampling:
% nyquistSampling=1 corresponds to 2 pixels per fwhm
la.nyquistSampling = 0.5;

%% increasing the field of view
% fieldStopSize is given in units of fwhm, for example for a nxn input wave
% and nyquistSampling=1, the default value of fieldStopSize is
% (n/nLenslet)/2
la.fieldStopSize = 6;

%% back to original sampling and field of view
la.nyquistSampling = 1;
la.fieldStopSize = 3;

%% a random aberration
zern = zernike(5:6,'resolution',nPx);
zern.c = 4*rand(2,1);
zern.lex = false;
la.wave = zern.wave;
propagateThrough(la)

%% a random aberration function
f = @(x) zern.wave;
la.wave = f;
for k=1:100
    o = (k-1).*2*pi/99;
    zern.c = 4.*[cos(o);sin(o)];
    propagateThrough(la)
    drawnow
end

%% stacked waves
la.wave = cat(3,utilities.piston(nPx),zern.wave);
la.imageletsListener.Enabled = false;
propagateThrough(la)
figure
imagesc( [la.imagelets(:,:,1) , la.imagelets(:,:,2)] )
axis equal tight
colorbar('location','NorthOutside')

%% sum stacked waves
zern = zernike(4,'resolution',nPx);
zern.lex = false;
zern.c = 6*rand;
wave1  = zern.wave;
zern.c = -6*rand;
wave2  = zern.wave;
la.wave = cat(3,wave1,wave2);
la.sumStack = true;
propagateThrough(la)
imagesc( la )
la.sumStack = false;

%% <matlab:doc('lensletArray') lensletArray> combined with an <matlab:doc('telescope') telescope> object
atm = atmosphere(2.2e-6,1.5,25,...
    'altitude',2e3,...
    'fractionnalR0',1,...
    'windSpeed',10,...
    'windDirection',0);
sys = telescope(30,...
    'fieldOfViewInArcMin',2,...
    'resolution',nPx*3,...
    'samplingTime',1/100,...
    'opticalAberration',atm);
update(sys)
imagesc(sys)
la.lightSource = source; % NGS
la.wave = sys;
propagateThrough(la)
imagesc( la )
la.imageletsListener.Enabled = true;
for k=1:50
    propagateThrough(la)
    drawnow
end
%%
% a simple Laser Guide Star
la.lightSource = source('height',[89.6,90,90.4].*1e3,'wavelength',589e-9); 
la.sumStack = true;
sys.focalDistance = 90e3;
propagateThrough(la)
for k=1:50
    propagateThrough(la)
    drawnow
end

