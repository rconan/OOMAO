%% UTILITIES HOWTO
% Demonstrate the use of the  <matlab:doc('utilities') utilities> class
%%

%% piston
% Create a circular or square mask with 1 inside the maks and 0 outside
%%
% * a disc in a 32x32 array
imagesc(utilities.piston(32))
axis square
snapnow
%%
% * a centered 24 pixel diameter disc in a 32x32 array
imagesc(utilities.piston(24,32))
axis square
snapnow
%%
% * a offseted 24 pixel diameter disc in a 32x32 array
imagesc(utilities.piston(18,32,-4,4))
axis square
snapnow
%%
% * a offseted 24 pixel size square in a 32x32 array
imagesc(utilities.piston(18,32,-4,4,'shape','square'))
axis square
snapnow

%% toggleFrame
% Reshape an array from 2 dimensions to 3 dmensions and vice-versa
%%
% a 3d array
a = cat( 3 , magic(4) , round(magic(4)/2) )
%%
% reshape into a 2d array
a2d = utilities.toggleFrame(a)
%%
% and back to 3d
a3D = utilities.toggleFrame(a2d)
%%
% forcing 2d reshape
utilities.toggleFrame(a,2)
%%
% forcing 3d reshape
utilities.toggleFrame(a,3)

%% rearrange
% Creates array index to reshape array as follows
a = reshape(1:36,6,6)
index = utilities.rearrange(size(a),[3,3])
utilities.rearrange(size(a),[3,3],[],'column')
reshape( a(index) ,[3,3,4])