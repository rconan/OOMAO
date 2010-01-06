%% DETECTOR HOWTO
% Demonstrate the use of the <matlab:doc('detector') detector> class

%% A camera for a Shach-Hartmann wavefront sensor

%%
% the lenslet array
nLenslet = 20;
nPxPerLenset = 5;
la = lensletArray(nLenslet);

%%
% the telescope
res = nLenslet*nPxPerLenset;
cfht = telescope(3.6,'resolution',res);

%%
% the camera
cam = detector(res)
cam.frameGrabber = la;
source*cfht*la;
grab(cam)
imagesc(cam)
cam.frameListener.Enabled = true;

%%
% adding noise
cam.readOutNoise = 0.1;
grab(cam)


%%
% camera free run mode
cam.paceMaker.Period = 0.5;
start(cam.paceMaker)
pause(5)
stop(cam.paceMaker)


