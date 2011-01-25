function varargout = bilinearSplineInterpMat(gs,atm,tel,groundGrid)
%% BILINEARSPLINEINTERP Bilinear interpolation matrix
%
% [H,mask] = bilinearSplineInterpMat(gs,atm,tel,groundGrid) computes the
% bilinear interpolation matrices for a given system made of source,
% atmosphere and telescope objects and from the pupil sampling grid at the
% ground

nLayer      = atm.nLayer;
nGs         = length(gs);
nPxGround   = length(groundGrid);

% Ground sampling
% [xGround,yGround] ...
%            = utilities.cartAndPol(nPxGround,tel.R);
[xGround,yGround] ...
           = meshgrid(linspace(-1,1,nPxGround)*tel.R);
xGround    = xGround(groundGrid);
yGround    = yGround(groundGrid);
groundGrid = groundGrid(:)';

pitchGround = tel.D/(nPxGround-1);

% cell of bi-linear interpolation operators
H    = cell(nGs , nLayer);
% cell of phase layer mask where the bilinear spline are non-zeros
mask = cell(1   , nLayer);

fprintf('___ BI-HARMONIC OPERATOR ___\n')

for kLayer = 1:nLayer
    
    % Layer sampling
    D = atm.layer(kLayer).D;
%     nPxLayer        = ceil(D/pitchGround) + 1;
    nPxLayer        = floor(D/pitchGround) + 1;
    newD = (nPxLayer-1)*pitchGround;
    while newD<D
        nPxLayer        = nPxLayer + 2;
        newD = (nPxLayer-1)*pitchGround;
    end       
    D = newD;
    [xLayer,yLayer] = utilities.cartAndPol(nPxLayer,D/2);
%     [xLayer,yLayer] = meshgrid(linspace(-1,1,nPxLayer)*D/2);

    mask{kLayer} = false;
    
    for kGs=1:nGs
        
        fprintf(' [%d,%d]',kGs,kLayer)
        
        height = atm.layer(kLayer).altitude;
        % pupil center in layer
        beta   = gs(kGs).directionVector*height;
        
        if height==0
            mask{kLayer} = mask{kLayer} | groundGrid;
            nH = numel(xGround);
            mH = nPxGround^2;
            iH = 1:nH;
            jH = find(groundGrid);
            sH = ones(1,length(jH));
            H{kGs,kLayer} = sparse(iH,jH,sH,nH,mH);
        else
            H{kGs,kLayer} = bilinearSplineInterp(xLayer,yLayer,pitchGround,xGround - beta(1),yGround - beta(2));
            mask{kLayer} = mask{kLayer} | ( ~all(H{kGs,kLayer}==0) );
        end
        
    end
    
    fprintf('\n')
    
end

fprintf('----------------------------\n')

varargout{1} = H;
varargout{2} = mask;

end

function H = bilinearSplineInterp(xo,yo,do,xi,yi)
%% BILINEARSPLINEINTERP Bilinear interpolation
%
% H = bilinearSplineInterp(xo,yo,do,xi,yi,d); computes the sparse matrix H
% to perform the bilinear interpolation zi = H*zo, where zo and zi are
% defined on the meshes [xo;yo] and [xi;yi], respectively. do is the mesh
% step size of [xo;yo]

xo = xo(:)';
yo = yo(:)';
xi = xi(:);
yi = yi(:);

ni = length(xi);
if ni~=length(yi)
    error('xi and yi must have the same length.')
end
no = length(xo);
if no~=length(yo)
    error('xo and yo must have the same length.')
end

u = bsxfun(@minus,xi,xo)/do;
v = bsxfun(@minus,yi,yo)/do;

H = linearSpline(u).*linearSpline(v);

% [i,j,s] = find(H);
% [m,n]   = size(H);
% as      = abs(s);
% index   = as > eps(max(as)); % round-off error filter
%
% H = sparse(i(index),j(index),s(index),m,n);


    function y = linearSpline(x)
        %% LINEARSPLINE Linear spline function
        %
        % y = linearSpline(x) computes the function y = 1 - |x| for |x|<1 and y = 0
        % elsewhere
        
        [m,n] = size(x);
        x     = abs(x);
        index = x < 1;
        [i,j] = find(index);
        s     = 1 - x(index);
        y = sparse(i,j,s,m,n);
        
    end

end
