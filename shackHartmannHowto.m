%% SHACKHARTMANN HOWTO 
% Demonstrate the use of the <matlab:doc('shackHartmann') shackHartmann> class

%% shackHartmann definition
wfs = shackHartmann(20,120,0.75);
setValidLenslet(wfs,utilities.piston(120))
tel = telescope(8,'resolution',120);
ngs = source;
wfs.lenslets.lightSource = ngs;
wfs.lenslets.wave = tel;
grabAndProcess(wfs)
wfs.referenceSlopes = wfs.slopes;
dataProcessing(wfs)
figure
slopesAndFrameDisplay(wfs)
wfs.camera.frameListener.Enabled = true;
wfs.slopesListener.Enabled = true;
% dataProcessing(wfs)

%% a random aberration
zern = zernike(5:6,'resolution',120,'pupil',tel.pupil);
zern.c = 20*(2*rand(zern.nMode,1)-1);
zern.lex = false;
wfs.lenslets.wave = zern.wave;
grabAndProcess(wfs)

%% a random aberration function
f = @(x) zern.wave;
wfs.lenslets.wave = f;
for k=1:100
    o = (k-1).*2*pi/99;
    zern.c = 10.*[cos(o);sin(o)];
    grabAndProcess(wfs)
    drawnow
end

