clear
ngs = source;
nPx = 60;
nLenslet = 10;
tel = telescope(1,'resolution',nPx);

zern = zernike(tel,2);
zern.c = ngs.wavelength/10;
pyr = pyramid(nLenslet,nPx,'modulation',4,'binning',1,'c',2);
%pyr.modulation = 1;
%pyr.binning = 2;

ngs.*tel*pyr;
%pyr.INIT

%%
ngs.*tel*pyr;

figure(103)
subplot(2,1,1)
imagesc(pyr.camera)
axis equal tight
subplot(2,1,2)
slopesDisplay(pyr)
axis equal tight

pyr.camera.frameListener = true;
pyr.slopesListener = true;

%%
ast = source('asterism',{[3,arcsec(30),0]});
ast.*tel*zern*pyr;
imagesc(pyr.camera)
slopesDisplay(pyr)

%%
bif = influenceFunction('monotonic',0.5);
dm = deformableMirror(nLenslet+1,'modes',bif,'resolution',nPx,...
    'validActuator',pyr.validActuator);
ngs.*tel;
calib = calibration(dm,pyr,ngs,-ngs.wavelength/80);

%for k=1:dm.nValidActuator;dm.coefs(k)=ngs.wavelength/40;+ngs;drawnow;dm.coefs(k)=0;end
