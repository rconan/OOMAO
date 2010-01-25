%% SHACKHARTMANN HOWTO 
% Demonstrate the use of the <matlab:doc('shackHartmann') shackHartmann> class

%% shackHartmann definition
wfs = shackHartmann(20,120,0.75);
setValidLenslet(wfs,utilities.piston(120))
tel = telescope(8,'resolution',120);
ngs = source;
% wfs.lenslets.lightSource = ngs;
% relay(tel,ngs);
ngs = ngs.*tel*wfs.lenslets;
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
% reset(ngs)*tel*zern*wfs.lenslets;
% grabAndProcess(wfs)
ngs = ngs.*tel*zern*wfs;

%% a random aberration function
for k=1:100
    o = (k-1).*2*pi/99;
    zern.c = 10.*[cos(o);sin(o)];
    +ngs;
    drawnow
end

