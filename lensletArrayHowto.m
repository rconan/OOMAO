%% LENSLETARRAY HOWTO
% Demonstrate the use of the <matlab:doc('lensletArray') lensletArray> class

%% A simple 10x10 lenslet array
nlenslet = 10;
nPx = 60;
la = lensletArray(nlenslet)
% a circular pupil
ngs = source*utilities.piston(nPx);
% a circular pupil propagated throught the lenslet array
propagateThrough(la,ngs)
% and displayed
imagesc(la)
% auto-update of the display
la.imageletsListener.Enabled = true;

%% changing the sampling:
% nyquistSampling=1 corresponds to 2 pixels per fwhm
la.nyquistSampling = 0.5;
propagateThrough(la,ngs)

%% increasing the field of view
% fieldStopSize is given in units of fwhm, for example for a nXn input wave
% and nyquistSampling=1, the default value of fieldStopSize is
% (n/nLenslet)/2
la.fieldStopSize = 6;
propagateThrough(la,ngs)

%% back to original sampling and field of view
la.nyquistSampling = 1;
la.fieldStopSize = 3;
propagateThrough(la,ngs)

%% a random aberration
zern = zernike(5:6,'resolution',nPx);
zern.c = 4*rand(2,1);
propagateThrough(la,ngs*zern)

%% a random aberration function
for k=1:100
    o = (k-1).*2*pi/99;
    zern.c = 4.*[cos(o);sin(o)];
    propagateThrough(la,reset(ngs)*zern)
    drawnow
end

%% stacked waves
zern.lex = false;
ngs = source*cat(3,utilities.piston(nPx),zern.wave);
la.imageletsListener.Enabled = false;
propagateThrough(la,ngs)
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
ngs = source*cat(3,wave1,wave2);
la.sumStack = true;
propagateThrough(la,ngs)
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
%%
ngs = source*sys;
propagateThrough(la,ngs)
imagesc( la )
la.imageletsListener.Enabled = true;
for k=1:50
    reset(ngs)*update(sys)*la;
    drawnow
end
%%
% a simple Laser Guide Star
lgs = source('height',[89.6,90,90.4].*1e3,'wavelength',589e-9); 
sys.focalDistance = 90e3;
propagateThrough(la,lgs*sys)
for k=1:50
    reset(lgs)*update(sys)*la;
    drawnow
end

