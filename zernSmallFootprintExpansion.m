% Zernike small footprint expansion

Dratio = 3;
nradialOrder ...
       = 5;
nPx    = 32;
zern   = zernike(1:zernike.nModeFromRadialOrder(nradialOrder),nPx);
deltaX = 0.5;
deltaY = 0.5;
[x,y,r,o] = utilities.cartAndPol(nPx);

fun = @(r,o,zi,ni,mi,zj,nj,mj) ...
    r.*zernike.funOld(zi,ni,mi,...
    complex( (r.*cos(o)+deltaX)/Dratio,(r.*sin(o)+deltaY)/Dratio ) ) .* ...
    zernike.funOld(zj,nj,mj,...
    r.*exp(1i.*o)) .* (r<=1);
zernBigCut = @(r,o,zi,ni,mi) ...
    zernike.funOld(zi,ni,mi,...
    complex( (r.*cos(o)+deltaX)/Dratio,(r.*sin(o)+deltaY)/Dratio ) ).* (r<=1);

nMode = length(zern.j);
zi = 8;
ni = 3;
mi = 1;
znmi = repmat({[zi;ni;mi]},1,nMode);
znmj = mat2cell([zern.j;zern.n;zern.m],3,ones(1,nMode));
znmj = repmat(znmj,nMode,1);
znmi = znmj';
index = tril(true(nMode));
pij = @(znmi,znmj) quad2d(@(r,o) fun(r,o,znmi(1),znmi(2),znmi(3),znmj(1),znmj(2),znmj(3)) ,0,1,0,2*pi);
P = zeros(nMode);
tic
P(index) = cellfun(pij,znmi(index),znmj(index))./pi;
toc
            P(abs(P)<1e-6) = 0;
% zern.c = P';    
% zern0   = zernike(zi,nPx);
% zern.lex = false;
% zern0.lex = false;
% figure(1)
% subplot(2,2,1)
% imagesc(zern0.phase)
% axis square
% colorbar
% subplot(2,2,2)
% imagesc(zernBigCut(r,o,zi,ni,mi))
% axis square
% colorbar
% subplot(2,2,[3,4])
% imagesc([zernBigCut(r,o,zi,ni,mi),zern.phase])
% axis equal tight
% colorbar