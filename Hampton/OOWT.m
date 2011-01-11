% Much of the following code is copied from Rodolphe Conan's example
% "adaptiveOpticsHowTo.m". His comments have been removed here and can
% still be found there. Comments in this code show where my program differ
% from his.
%
% The challenge of this program is to show that I can program a woofer
% tweeter system with the tools developed by Dr. Conan.
%
% Peter Hampton, August 4, 2010

atm = atmosphere(photometry.V,0.15,30,...
    'altitude',[0,4,10]*1e3,...
    'fractionnalR0',[0.7,0.25,0.05],...
    'windSpeed',[5,10,20],...
    'windDirection',[0,pi/4,pi]);

nLenslet = 15;

nPx_Lenslet = 6;

sampleTime = 1/100;

% nPx is entirely dependant on the choice of the number of lenslets and the
% choice of the number of pixels per lenslet. Attempted nLenslet = 31 but
% that is prohibitively computationally expensive.

nPx = nPx_Lenslet*nLenslet;

tel = telescope(3.6,...
    'fieldOfViewInArcMin',2.5,...
    'resolution',nPx,...
    'samplingTime',1/100);

ngs = source('wavelength',photometry.J,'height',90000,'magnitude',[]);

tel.focalDistance = ngs.height;

wfs = shackHartmann(nLenslet,nPx,0.5);

ngs = ngs.*tel;
ngs = ngs*wfs;

setValidLenslet(wfs)

+wfs;

wfs.referenceSlopes = wfs.slopes;

+wfs;

tmbif = influenceFunction('monotonic',25/100);

wmbif = influenceFunction('monotonic',25/100);

nActuator = nLenslet + 1;

% Woofer (wm) is chosen to be 8 x 8. The tweeter is tm. These replace the
% single instance of dm.

wmvAct = ones(8);
wmvAct(1,[1 2 7 8]) = 0;
wmvAct(8,[1 2 7 8]) = 0;
wmvAct(2,[1 8]) = 0;
wmvAct(7,[1 8]) = 0;
wmvAct = logical(wmvAct);

wm = deformableMirror(8,...
    'modes',wmbif,...
    'resolution',nPx,...
    'validActuator',wmvAct); % 

tm = deformableMirror(nActuator,...
    'modes',tmbif,...
    'resolution',nPx,...
    'validActuator',wfs.validActuator,...
    'distortionTau', sampleTime,...
    'distortionSaturation',3,...
    'distortionConstant',2,...
    'sampleTime',sampleTime);

stroke = ngs.wavelength/2;

tm.coefs = eye(tm.nValidActuator)*stroke;

wm.coefs = eye(wm.nValidActuator)*stroke;



% Leave the tm flat and project each woofer actuator onto the wfs to obtain
% the woofer's interaction matrix

ngs=ngs.*tel;
ngs=ngs*wm;
ngs=ngs*wfs;

wmInteractionMatrix = wfs.slopes./stroke;

% Leave the wm flat and project each tweeter actuator onto the wfs to obtain
% the tweeter's interaction matrix

ngs=ngs.*tel;
ngs=ngs*tm;
ngs=ngs*wfs;

tmInteractionMatrix = wfs.slopes./stroke;

[U,S,V]=svd(tmInteractionMatrix);
s=diag(S);

iS = diag(1./s(find(s > 0.1*s(1))));
[nS,nC] = size(tmInteractionMatrix);
iS(nS,nC) = 0;

tmRecon = V*iS'*U';

[Utm,S,Vwm] = svd(tmRecon*wmInteractionMatrix);

Utm = Utm(:,1:wm.nValidActuator);
Sinv = diag(1./diag(S));

tel=tel+atm;
wm.coefs = 0;
tm.coefs = 0;

g = 0.5;
    ngs=ngs.*+tel;  
    ngs=ngs*wm*tm*wfs;
for kIteration = 1:200
    
    ngs=ngs.*+tel;          % arriving aberrated wave front

    turbPhase = ngs.meanRmPhase;
    
    dmCoefs = tmRecon*wfs.slopes;   % this line is before propogation to simulate the delay
    
    ngs=ngs*wm*tm*wfs;
%    total(kIteration) = var(ngs);

 

%    residue(kIteration) = var(ngs);

   
    % Integrating the DM coefficients
    
    wmModes = Utm'*dmCoefs;
    
    tmCoefs = dmCoefs - Utm*wmModes;
    
    wmCoefs = Vwm*Sinv*wmModes;
    
    tm.coefs = tm.coefs - g*tmCoefs;
    
    wm.coefs = wm.coefs - g*wmCoefs;
    
    
    % Display of turbulence and residual phase
    figure(1)
    title(kIteration)
    subplot(211)
    imagesc([turbPhase,ngs.meanRmPhase]);
    subplot(212)
    imagesc([wm.surface,tm.surface]);
    drawnow
    
    
    
end




figure(10)
subplot(121)
imagesc(wmInteractionMatrix)
xlabel('WM actuators')
ylabel('WFS slopes [px]')
ylabel(colorbar,'slopes/actuator stroke')

subplot(122)
imagesc(tmInteractionMatrix)
xlabel('TM actuators')
ylabel('WFS slopes [px]')
ylabel(colorbar,'slopes/actuator stroke')