%% OOMAO TUTORIAL
close all
clear

%%%
% NGS AO
%%%

%% telescope class (mandatory field: the diameter)
tel = telescope(1,'resolution',100,'fieldOfViewInArcsec',30,'samplingTime',1/500);
disp(tel)

%% wavefront sensor class (mandatory field: the size of the lenslet array
% (10x10 on the following), the camera resolution)
wfs = shackHartmann(10,100,0.85);

%% source class
ngs =source;

%% propagation of source
ngs = ngs.*tel*wfs;

%% wfs initialization: set the mask of valid lenslet and the reference
% slopes
wfs.INIT
+wfs
figure
imagesc(wfs.camera)

%% deformable mirror influence function class
bifa = influenceFunction('monotonic',0.75);
figure,show(bifa)
bifb = influenceFunction('overshoot',0.75);
figure,show(bifb)

%% deformable mirror class
dm = deformableMirror(11,'modes',bifa,'resolution',tel.resolution,'validActuator',wfs.validActuator);

%% dm calibration (method 1)
dm.coefs = eye(dm.nValidActuator)*ngs.wavelength;
ngs = ngs.*tel*dm*wfs;
figure,imagesc(wfs.slopes),colorbar

%% dm calibration (method 2)
ngs = ngs.*tel;
calibDm = calibration(dm,wfs,ngs,ngs.wavelength);
calibDm.threshold = 1e6;
disp(calibDm)

%% atmosphere class
atm = atmosphere(photometry.V,15e-2,30,'altitude',5e3,'windSpeed',10,'windDirection',pi/3);
tel = tel + atm;
figure
imagesc(tel)

%% phase screens Taylor (frozen flow) translation
while true,imagesc(+tel),drawnow,end

%% turbulence analysis with wavefront sensor
wfs.camera.frameListener.Enabled = true;
ngs = ngs.*tel*wfs;
while true,imagesc(+tel),+ngs,drawnow,end

%% photon noise on wfs detector
wfs.camera.photonNoise = true;
ngs.magnitude = 10;

%% readout-noise on wfs detector
ngs.magnitude = 6;
wfs.camera.readOutNoise = 5;

%% resetting the detector to a noiseless one
wfs.camera.photonNoise = false;
wfs.camera.readOutNoise = 0;
+wfs

%% imaging camera
cam = imager(100);
science = source('wavelength',photometry.K);
science = science.*tel*cam;
figure
imagesc(cam)
cam.frameListener.Enabled = true;
while true;+tel;+science;drawnow;end

%% initializing imaging camera for Strehl ratio computing
tel = tel - atm;
+science
cam.referenceFrame = cam.frame;
+science
cam.strehl
tel = tel + atm;
+science
cam.strehl

%% long exposure imaging camera
cam.exposureTime = 1000;
cam.clockRate = 1;
flush(cam)
cam.strehl
flush(cam)

%% closed-loop adaptive optics
gain = 0.5;
dm.coefs = 0;
ngs.logging = true;
ngs = ngs.*tel*dm*wfs;
figure
h = imagesc(ngs.meanRmOpd*1e6);
colorbar
cam.exposureTime = 150;
cam.clockRate = 1;
science = science.*tel*dm*cam;
dmCoefs = size(dm.coefs,2);
pause(1)
for k=1:150
    +tel
    +ngs
    +science
    dm.coefs = dm.coefs - gain*calibDm.M*wfs.slopes;
%     dmCoefs(:,2) = dmCoefs(:,1) - gain*calibDm.M*wfs.slopes;
%     dm.coefs = dmCoefs(:,1);
%     dmCoefs(:,1) = dmCoefs(:,2);
    set(h,'Cdata',ngs.meanRmOpd*1e6)
    drawnow
end

%% reporting performance
hf = figure;
plot(1e6*sqrt(ngs.phaseVar(1:150*2))/ngs.waveNumber,'.')

%% open-loop adaptive optics
dm.coefs = 0;
ngs = ngs.*tel*wfs*dm;
figure
h = imagesc(ngs.meanRmOpd*1e6);
colorbar
pause(1)
for k=1:150
    +tel
    +ngs
    dm.coefs = -calibDm.M*wfs.slopes;
    set(h,'Cdata',ngs.meanRmOpd*1e6)
    drawnow
end

%% reporting performance
figure(hf)
hold all
plot(1e6*sqrt(ngs.phaseVar((1:150*2)+300))/ngs.waveNumber,'.')
results = reshape(ngs.phaseVar(1:600),300,2);
fprintf('mean residual wfe in closed-loop: %4.2fnm and in open-loop: %4.2fnm\n',1e9*sqrt(mean(results(20:end,:)))/ngs.waveNumber)

%%%
% MULTIPLE SOURCE
%%%

%% binary system
binary = source('zenith',arcsec([0,15]),'azimuth',[0,0]);
figure,imagesc(tel,binary)
binary = binary.*tel*wfs;

%% laser guide star (point source)
lgs = source('asterism',{[3,arcsec(15),0]},'height',90e3);
figure,imagesc(tel,lgs)
lgs = lgs.*tel*wfs;
tel = tel - atm;
+lgs

%% laser guide star (point source defocus)
set(lgs,'height',92e3)
+lgs

%% laser guide star (point source defocus and off-axis launch)
set(lgs,'viewPoint',[25,0])
+lgs

%% extended laser guide star (point source)
lgs = source('asterism',{[3,arcsec(15),0]},'height',[88,92]*1e3);
set(lgs,'objectiveFocalLength',90e3)
lgs = lgs.*tel*wfs;

%% extended laser guide star (off-axis launch)
set(lgs,'viewPoint',[25,0])
+lgs

%% extended laser guide star (3 off-axis launch)
% set(lgs(1,1,:),'viewPoint',[0,0])
% set(lgs(1,2,:),'viewPoint',[10,10])
% set(lgs(1,3,:),'viewPoint',[0,25])
% +lgs
% 
% %% extended laser guide star (3 off-axis launch and z-profile)
% lgs(1,2,1).nPhoton = 0;
% +lgs
% %%
% lgs(1,2,1).nPhoton = lgs(1,2,2).nPhoton*2;
% +lgs
% %%
% lgs(1,2,1).nPhoton = lgs(1,2,2).nPhoton/2;
% +lgs

