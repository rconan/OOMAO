function varargout = sparseInterpMatrix(inGridMask,outGridMask,tel,src,atm)

nSrc    = size(src,2);
H       = cell(nSrc,atm.nLayer);
inGrid  = cell(nSrc,atm.nLayer);
outGrid = cell(nSrc,atm.nLayer);
interpPhase ...
        = cell(nSrc,atm.nLayer);
if ~iscell(inGridMask)
    inGridMask = repmat( {inGridMask} , 1 , atm.nLayer);
end

npMap = length(outGridMask);
indOutMap = find(outGridMask);
numelOutMap = length(indOutMap);

fprintf('___ Sparse Interpolation Matrix (%dX%d) ___\n',nSrc,atm.nLayer) 
for kSrc=1:nSrc
    for kLayer = 1:atm.nLayer
                
        fprintf(' [%2d,%2d]',kSrc,kLayer)
        
        nMap  = length(inGridMask{kLayer});
        indInMap = find(inGridMask{kLayer});
        numelInMap = length(indInMap);
        % First 2x2 block
        ind0 = sub2ind([nMap,nMap],[1 2 1 2],[1 1 2 2])';
        % First column of 2x2 block
        nBlock = max(nMap - 1,2);
        ind = bsxfun(@plus,ind0,0:nMap-2);
        % Other columns of 2x2 block
        if nMap>2
            ind = bsxfun(@plus,ind(:),0:nMap:nMap*(nBlock-1));
            ind = reshape(ind,4,[]);
        end

        % interp. grid
        alpha = (1-atm.layer(kLayer).altitude/src(kSrc).height).*...
            tel.D/atm.layer(kLayer).D;
        beta = atm.layer(kLayer).altitude.*[ src(kSrc).directionVector.x , src(kSrc).directionVector.y ].*...
            nMap./atm.layer(kLayer).D;
        [x,y]   = meshgrid((1-nMap:2:nMap-1)/2);
        [xp,yp] = meshgrid(((nMap-1)/(npMap-1))*(1-npMap:2:npMap-1)/2);
        xp = alpha*xp + beta(1);
        yp = alpha*yp + beta(2);
        inGrid{kSrc,kLayer}  = {x,y};
        outGrid{kSrc,kLayer} = {xp,yp};
%         interpPhase{kSrc,kLayer} = ...
%             griddata(x,y,atm.layer(kLayer).phase,xp,yp);
        % origin shift
        xp = xp - x(1);
        x  = x  - x(1);
        yp = yp - y(1);
        y  = y  - y(1);

        xx = x(ind);
        yy = y(ind);
        
        u   = 1:4;
        j_x = zeros(1,4*numelOutMap);
        w   = zeros(1,4*numelOutMap);

        for kMap=indOutMap'
            
            % sum of distance squared from 2x2 block vertices to point (xp,yp)
            rho2 = sum((xx - xp(kMap)).^2 + (yy - yp(kMap)).^2);
            % the 2x2 block
            indp = ind(:,rho2==min(rho2));
            % interpolation
            j_x(u) = indp(:,1);
            %     fprintf('(xp,yp)=(%3.1f,%3.1f), (x0,y0)=(%3.1f,%3.1f)\n',...
            %         xp(kMap),yp(kMap),x(indp(1,1)),y(indp(1,1)))
            % origin shift
            xi = xp(kMap) - x(indp(1,1));
            yi = yp(kMap) - y(indp(1,1));
            % interp. weight
            w(u(1)) = (1-xi)*(1-yi);
            w(u(2)) = xi*(1-yi);
            w(u(4)) = xi*yi;
            w(u(3)) = (1-xi)*yi;
            %     fprintf('(xi,yi,zi)=(%3.1f,%3.1f,%3.1f), w=(%3.1f,%3.1f,%3.1f,%3.1f), z=(%3.1f,%3.1f,%3.1f,%3.1f)\n',...
            %         xi,yi,sum(w(u).*z(j_x(u))),w(u),z(j_x(u)))
            u = u + 4;
            
        end
        i_x = 1:numelOutMap;
        i_x = i_x(ones(1,4),:);
        H{kSrc,kLayer} = sparse(i_x,j_x,w,numelOutMap,nMap^2);
        H{kSrc,kLayer}(:,~inGridMask{kLayer}) = [];

    end
    fprintf('\n')
end
fprintf('-----------------------------------\n')

varargout{1} = H;
if nargout>1
    varargout{2} = inGrid;
    varargout{3} = outGrid;
    varargout{4} = interpPhase;
end