% Reflect.m
%   Extends rectangular data to fill a square by reflecting the gradient
%   data
%   Peter Hampton
%   Copyright August 2008

function [dzdx dzdy] = Reflect(dzdx,dzdy)
[row col] = size(dzdx);
if row > 2047
    dzdx = dzdx(1:2047,:);
    dzdy = dzdy(1:2047,:);
end
if col > 2047
    dzdx = dzdx(:,1:2047);
    dzdy = dzdy(:,1:2047);
end
[row col] = size(dzdx);
M = ceil(log2(size(dzdx)+1));
dzdx(row+1:2^M(1),:) = dzdx(row:-1:1+2*row-2^M(1),:);
dzdx(:,col+1:2^M(2)) = -dzdx(:,col:-1:1+2*col-2^M(2));
dzdy(row+1:2^M(1),:) = -dzdy(row:-1:1+2*row-2^M(1),:);
dzdy(:,col+1:2^M(2)) = dzdy(:,col:-1:1+2*col-2^M(2));
%only one of the following for loops will execute
for k = M(1):M(2)-1
    dzdx(2^k+1:2^(k+1),:) = dzdx(2^k:-1:1,:);
    dzdy(2^k+1:2^(k+1),:) = -dzdy(2^k:-1:1,:);
end
for k = M(2):M(1)-1
    dzdx(:,2^k+1:2^(k+1)) = -dzdx(:,2^k:-1:1);
    dzdy(:,2^k+1:2^(k+1)) = dzdy(:,2^k:-1:1);
end