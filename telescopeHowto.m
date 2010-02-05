%% TELESCOPE HOWTO
% Demonstrate the use of the <matlab:doc('telescope') telescope> class
%%

%% 
% $\lambda=2.2$
% micron,
%  $r_0=80$
% cm, 
%  ${\cal L}_0=25$
% m and diameter(D)=10m
%% single layer case
% the layer is set a 2km, the wind blows a 10m/s in
% the X direction, the field of view is 2', the 10m pupil is sampled with
% 60pixels and the phase screens are sampled at 500Hz
% The <matlab:doc('AT.update') update> method move the phase screen of one
% time step and <matlab:doc('AT.imagesc') imagesc> display the phase screen
atm = atmosphere(2.2e-6,0.8,25,...
    'altitude',2e3,...
    'fractionnalR0',1,...
    'windSpeed',10,...
    'windDirection',0)
tel = telescope(10,...
    'fieldOfViewInArcMin',2,...
    'resolution',60,...
    'samplingTime',1/500)
tel = tel + atm;
+tel;
imagesc(tel)
%% 
% At any time, the phase screen of each layer is stored in the phase
% property of the <matlab:doc('turbulenceLayer') layer> object 
tel.opticalAberration.layer
%%
% For a given direction in the sky, the geometric propagation of the phase
% through the turbulence layers into the telescope pupil are computed with
% the method <matlab:doc('AT.getPhaseScreen') getPhaseScreen>, the
% direction in the sky is given by the <matlab:doc('AT.source') source>
% object
%%
% * a natural guide star on axis
src = source.*tel;
figure
imagesc(src.meanRmPhase)
axis equal tight, colorbar('NorthOutside')
%%
% * a natural guide star off axis
src = source('zenith',30*cougarConstants.arcsec2radian,'azimuth',pi/4) .* tel;
imagesc(src.phase)
axis equal tight, colorbar('NorthOutside')
%%
% * an asterism of 4 stars conjugated at infinity
src =  source('asterism',{[0,0],[3,60*cougarConstants.arcsec2radian,30]}) .* tel;
imagesc(src.catPhase)
axis equal tight, colorbar('NorthOutside')
%%
% * a single source laser guide star
src = source('zenith',0,'azimuth',0,'height',91e3,'wavelength',589e-9);
tel.focalDistance = 90e3;
src = src .* tel;
imagesc(src.phase)
axis equal tight, colorbar('NorthOutside')

%% a 3 layer case
delete(imagesc(tel))
atm = atmosphere(2.2e-6,0.8,25,...
    'altitude',[0,10,15].*1e3,...
    'fractionnalR0',[0.7,0.2,0.1],...
    'windSpeed',[10,5,15],...
    'windDirection',[0,pi/4,pi/2]);
tel=tel+atm;
imagesc(+tel)
%%
tel.opticalAberration.layer(1)
tel.opticalAberration.layer(2)
tel.opticalAberration.layer(3)
%%
% * a natural guide star on axis
src = source .* tel;
imagesc(src.phase)
axis equal tight, colorbar('NorthOutside')
%%
% * a natural guide star off axis
src = source('zenith',30*cougarConstants.arcsec2radian,'azimuth',pi/4) .* tel;
imagesc(src.phase)
axis equal tight, colorbar('NorthOutside')
%%
% * an asterism of 4 stars conjugated at infinity
src =  source('asterism',{[0,0],[3,60*cougarConstants.arcsec2radian,30]}) .* tel;
imagesc([src.meanRmPhase])
axis equal tight, colorbar('NorthOutside')
%%
% * a single source laser guide star
src = source('zenith',0,'azimuth',0,'height',91e3,'wavelength',589e-9);
tel.focalDistance = 90e3;
src =  src .* tel;
imagesc(src.phase)
axis equal tight, colorbar('NorthOutside')

