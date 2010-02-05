%% LENSLETARRAY HOWTO
% Demonstrate the use of the <matlab:doc('lensletArray') lensletArray> class

%% A simple 10x10 lenslet array
nlenslet = 10;
nPx = 60;
la = lensletArray(nlenslet)
% a circular pupil
ngs = source.*utilities.piston(nPx);
% a circular pupil propagated throught the lenslet array
ngs = ngs*la;
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
ngs = ngs.*zern*la;

%% a random aberration function
for k=1:100
    o = (k-1).*2*pi/99;
    zern.c = 4.*[cos(o);sin(o)];
    +ngs;
    drawnow
end

%% stacked waves
zern.lex = false;
ngs = source.*cat(3,utilities.piston(nPx),ngs.wave);
la.imageletsListener.Enabled = false;
ngs = ngs*la;
imagesc( [la.imagelets(:,:,1) , la.imagelets(:,:,2)] )
axis equal tight
colorbar('location','NorthOutside')

%% sum stacked waves
zern = zernike(4,'resolution',nPx);
zern.lex = false;
zern.c = 6*rand;
ngs1  = source.*zern;
zern.c = -6*rand;
ngs2  = source.*zern;
la.sumStack = true;
ngs = source.*cat(3,ngs1.wave,ngs2.wave)*la;
imagesc( la )
la.sumStack = false;

%% <matlab:doc('lensletArray') lensletArray> combined with an <matlab:doc('telescope') telescope> object
atm = atmosphere(2.2e-6,1.5,25,...
    'altitude',2e3,...
    'fractionnalR0',1,...
    'windSpeed',10,...
    'windDirection',0);
tel = telescope(30,...
    'fieldOfViewInArcMin',2,...
    'resolution',nPx*3,...
    'samplingTime',1/100)
tel = tel+atm;
imagesc(+tel)
%%
ngs = source.*tel*la;
imagesc( la )
la.imageletsListener.Enabled = true;
for k=1:50
    +tel;
    +ngs;
    drawnow
end
%%
% a simple Laser Guide Star
lgs = source('height',[89.6,90,90.4].*1e3,'wavelength',2.2e-6); 
tel.focalDistance = 90e3;
lgs = lgs.*tel*la;
for k=1:50
    +tel;
    +lgs;
    drawnow
end

