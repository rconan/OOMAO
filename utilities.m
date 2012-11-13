classdef utilities
    % UTILITIES Collection of various useful functions
    
    methods (Static)
        
        function out = piston(Npx,varargin)
            %% PISTON piston mode
            %
            % out = piston(Npx) Computes a piston on Npx pixel across the
            % diameter
            %
            % out = piston(Npx,nOut) Computes a piston on Npx pixel across
            % the diameter inside a square array of nOutXnOut pixels.
            %
            % out = piston(Npx,nOut,xOffset,yOffset) Computes a piston on
            % Npx pixel across the diameter inside a square array of
            % nOutXnOut pixels at xOffset and yOffset pixels from the
            % center.
            %
            % out = piston( ... ,'shape','square') By default the piston is
            % a disc but here it is forced to be a square
            %
            % out = piston( ... ,'shape','hex') By default the piston is
            % a disc but here it is forced to be hexagonal, nOut is equal
            % to twice the hexagonal side
            %
            % out = piston( ... ,'type','logical') By default the piston
            % values are in double but they can be casted into any types
            % supported by Matlab like logical
            
            p = inputParser;
            p.addRequired('Npx',@isnumeric);
            p.addOptional('nOut',Npx,@isnumeric);
            p.addOptional('xOffset',0,@isnumeric);
            p.addOptional('yOffset',0,@isnumeric);
            p.addParamValue('shape','disc',@ischar);
            p.addParamValue('type','double',@ischar);
            parse(p,Npx,varargin{:});
            param = p.Results;
            
            x = -(param.nOut-1)/2:(param.nOut-1)/2;
            u = x - param.xOffset;
            v = x - param.yOffset;
            
            [x,y,r,o] = utilities.cartAndPol(2.*u./Npx,2.*v./Npx);
            
            switch param.shape
                case 'disc'
                    out = double(r <= 1);
                case 'square'
                    out = double( abs(x)<=1 & abs(y)<=1 );
                case {'hex','hexagon'}
                    out = double( abs(x)<=sqrt(3)/2 & abs(y)<=x/sqrt(3)+1 & abs(y)<=-x/sqrt(3)+1 );
                otherwise
                    error('The piston shape is either a disc, a square or a hexagon')
            end
            
            switch param.type
                case 'logical'
                    out = logical(out);
                case 'double'
                otherwise
                    error('The piston type is either a double or logical')
            end
        end
        
        function varargout = cartAndPol(u,varargin)
            %% CARTANDPOL Cartesian and polar coordinate arrays
            %
            % [x,y,r,o] = cartAndPol(n) Shortcuts to u = ((1-n):2:(n-1))/n;
            % [x,y] = meshgrid(u);[o,r] = cart2pol(x,y);
            %
            % [x,y,r,o] = cartAndPol(n,R) Shortcuts to u = R*((1-n):2:(n-1))/n;
            % [x,y] = meshgrid(u);[o,r] = cart2pol(x,y);
            %
            % [x,y,r,o] = cartAndPol(u) Shortcuts to [x,y] = meshgrid(u);[o,r] =
            % cart2pol(x,y);
            %
            % [x,y,r,o] = cartAndPol(u,v) Shortcuts to [x,y] = meshgrid(u,v);[o,r] =
            % cart2pol(x,y);
            %
            % [x,y,r,o] = cartAndPol(n,[],type)
            %
            % [x,y,r,o] = cartAndPol(n,R,type)
            %
            % [x,y,r,o] = cartAndPol(u,[],type
            %
            % [x,y,r,o] = cartAndPol(u,v, type) Same as above but now the type of x, y,
            % r and o is now specified: double (default) or single
            %
            % [...] = cartAndPol(...,'offset',[xOffset,yOffset]) offsets
            % the grid by xOffset and yOffset
            %
            % [r,o] = cartAndPol(...,'output','polar')
            %
            % [r] = cartAndPol(...,'output','radius')
            
            
            p = inputParser;
            p.addRequired('u',@isnumeric);
            p.addOptional('v',[],@isnumeric);
            p.addParamValue('offset',[0,0],@isnumeric);
            p.addParamValue('type','double',@ischar);
            p.addParamValue('output','all',@ischar);
            p.parse(u, varargin{:})
            u      = p.Results.u;
            v      = p.Results.v;
            offset = p.Results.offset;
            type   = p.Results.type;
            output = p.Results.output;
            
            if isempty(v)
                if numel(u)==1
                    u = linspace(-1,1,u);%2*( -(u-1)/2:(u-1)/2 )/u;
                end
                v=u;
            elseif (numel(u)==1) && (numel(v)==1)
                u = linspace(-v,v,u);%2*v*( -(u-1)/2:(u-1)/2 )/u;%linspace(-v,v,u);
                v = u;
            end
            
            if strcmp(type,'single')
                [x,y] = meshgrid(single(u-offset(1)),single(v-offset(2)));
            else
                [x,y] = meshgrid(u-offset(1),v-offset(2));
            end
            [o,r] = cart2pol(x,y);
            
            switch output
                case 'all'
                    varargout{1} = x;
                    varargout{2} = y;
                    varargout{3} = r;
                    varargout{4} = o;
                case 'polar'
                    varargout{1} = r;
                    varargout{2} = o;
                case 'radius'
                    varargout{1} = r;
                otherwise
                    error('oomao:utilities:cartAndPol:wrongOutput',...
                        'Valid outputs are all, polar or radius.')
            end
            
        end
        
        function frame = toggleFrame(frame,toggle)
            %% TOGGLEFRAME 2D to 3D array reshaping
            %
            % out = toggleFrame(frame) reshapes a 2D array into a 3D array or a
            % 3d array into a 2D array
            %
            % out = toggleFrame(frame,toggle) reshapes the array into a 2D
            % array if toggle is equal to 2 or into a 3D array if toggle is
            % equal to 3
            
            n    = ndims(frame);
            dims = size(frame);
            if length(dims)==2
                dims(3) = 1;
            end
            
            if nargin<2
                if n==2
                    toggle = 3;
                else
                    toggle = 2;
                end
            end
            
            
            
            if n~=toggle || toggle==2
                switch toggle
                    case 2
                        %                         fprintf(' @(toggleFrame)> 2D: [%d,%d] !\n',dims(1)*dims(2),dims(3))
                        frame = reshape(frame,dims(1)*dims(2),dims(3));
                    case 3
                        m = sqrt(dims(1));
                        %                         fprintf(' @(toggleFrame)> 3D: [%d,%d,%d] !\n',m,m,dims(2))
                        frame = reshape(frame,[m,m,dims(2)]);
                end
            end
            
        end
        
        
        function index = rearrange(sizeArray,sizeSubArray,overlap,columnMajor)
            %5 REARRANGE Array linear index scrambling
            %
            % index = rearrange(sizeArray,sizeSubArray) Rearrange the linear index of
            % an array of size sizeArray in a 2D matrix where each row contains the
            % index of a sub-array of size sizeSubArray taken from the initial array
            %
            % index = rearrange(sizeArray,sizeSubArray,overlap) Rearrange the linear
            % index of  an array of size sizeArray in a 2D matrix where each row
            % contains the index of a sub-array of size sizeSubArray taken from the
            % initial array. The sub-arrays overlap a number overlap(1) of rows and
            % overlap(2) of columns.
            %
            % index = rearrange(sizeArray,sizeSubArray,[],'column') Same as above but
            % now the sub-array browse the array along the columns not the rows as
            % before.
            %
            % index = rearrange(sizeArray,sizeSubArray,overlap,'column') Same as above
            % with overlapping
            %
            % Example:
            % >> a = reshape(1:36,6,6)
            % a =
            %      1     7    13    19    25    31
            %      2     8    14    20    26    32
            %      3     9    15    21    27    33
            %      4    10    16    22    28    34
            %      5    11    17    23    29    35
            %      6    12    18    24    30    36
            % >> index = rearrange(size(a),[3,3])
            % ans =
            %      1     4    19    22
            %      2     5    20    23
            %      3     6    21    24
            %      7    10    25    28
            %      8    11    26    29
            %      9    12    27    30
            %     13    16    31    34
            %     14    17    32    35
            %     15    18    33    36
            % >> rearrange(size(a),[3,3],[],'column')
            % ans =
            %      1    19     4    22
            %      2    20     5    23
            %      3    21     6    24
            %      7    25    10    28
            %      8    26    11    29
            %      9    27    12    30
            %     13    31    16    34
            %     14    32    17    35
            %     15    33    18    36
            % >> reshape( a(index) ,[3,3,4])
            % ans(:,:,1) =
            %      1     7    13
            %      2     8    14
            %      3     9    15
            % ans(:,:,2) =
            %      4    10    16
            %      5    11    17
            %      6    12    18
            % ans(:,:,3) =
            %
            %     19    25    31
            %     20    26    32
            %     21    27    33
            % ans(:,:,4) =
            %     22    28    34
            %     23    29    35
            %     24    30    36
            
            % $Id: rearrange.m 409 2006-07-12 16:49:24Z aoteam $
            
            if nargin<3 || isempty(overlap)
                overlap = zeros(1,2);
            end
            
            n     = sizeArray(1);
            m     = sizeArray(2);
            k     = prod(sizeArray(3:end));
            nSub  = sizeSubArray(1);
            if numel(sizeSubArray)==1
                mSub = nSub;
            else
                mSub  = sizeSubArray(2);
            end
            
            if rem(n,2)
                % Odd n
                nNSub = (n+overlap(1))/nSub;
                mMSub = (m+overlap(2))/mSub;
            else
                % Even n
                nNSub = n/nSub + overlap(1);
                mMSub = m/mSub + overlap(2);
            end
            
            % Type the index array as an unsigned integer with the coding depending on
            % the value of the largest elements of the index array
            switch find(2.^(2.^(3:6))-1 > prod(sizeArray),1)
                case 1
                    uint = @(x) uint8(x);
                case 2
                    uint = @(x) uint16(x);
                case 3
                    uint = @(x) uint32(x);
                case 4
                    uint = @(x) uint64(x);
                otherwise
                    error('Array size to big')
            end
            
            % Sub-array index in array
            [i,j] = ndgrid(uint(1:nSub),uint(1:mSub));
            index = repmat( sub2ind( [n,m] , i(:) , j(:) ) , [ 1 , nNSub*mMSub*k ] );
            
            % Step index
            indexStep = ...
                repmat( uint(0:nNSub-1).'*(nSub-overlap(1))   , [  1   , mMSub*k ] ) + ...
                repmat( uint(0:mMSub*k-1)*(mSub-overlap(2))*n , [nNSub ,    1    ] );
            
            if nargin==4
                % Column major propagation of sub-array
                indexStep = indexStep.';
            end
            
            indexStep = repmat( reshape( indexStep , [1,nNSub*mMSub*k] ) , [nSub*mSub,1] );
            
            index = index + indexStep;
        end
        
        function out = sombrero(n,x)
            %% SOMBRERO Order n sombrero function
            %
            % out = sombrero(n,x) computes besselj(n,x)/x
            
            if n==0
                out = besselj(0,x)./x;
            else
                if n>1
                    out = zeros(size(x));
                else
                    out = 0.5*ones(size(x));
                end
                u = x~=0;
                x = x(u);
                out(u) = besselj(n,x)./x;
            end
        end
        
        function out = sinc(x)
            %% SINC Sinus cardinal function
            %
            % out = sinc(x) computes sin(pi*x)/(pi*x)
           
            out = ones(size(x));
            u = x~=0;
            x = x(u);
            out(u) = sin(pi*x)./(pi*x);
        end
        
        function out = fittingError(tel,atm,dm)
            %% FITTINGERROR Deformable mirror fitting error variance
            %
            % out = fittingError(telAtm,dm) computes the fitting error
            % variance of a a deformableMirror object for given telescope
            % and atmosphere objects
            
            c = (3/5)*(gamma(11/6)^2/pi^(8/3))*(24*gamma(6/5)/5)^(5/6);
            out = c*(tel.D/atm.r0)^(5/3)*...
                (dm.nValidActuator/pi + (tel.D/atm.L0)^2)^(-5/6);
        end
        
        function out = binning(frame,outRes)
            %% BINNING Frame binning
            % 
            % out = binning(frame,[n,m]) bins the frame pixels into a nXm
            % array; frame can be either a single frame or a data cube
            
            [n,m,nFrame] = size(frame);
            out          = zeros(outRes(1),outRes(2),nFrame);
            n1 = n/outRes(1);
            m2 = m/outRes(2);
            if n1==1 && m2==1
                out = frame;
                return
            end
            if n1==1
                for kFrame=1:nFrame
                    out(:,:,kFrame) = ...
                        reshape( ...
                        sum( ...
                        reshape( ...
                        frame(:,:,kFrame).', m2 , [] ) ...
                        ).' , ...
                        outRes(2) , [] ).';
                end
            elseif m2==1
                for kFrame=1:nFrame
                    out(:,:,kFrame) = ...
                        reshape( ...
                        sum( ...
                        reshape( frame(:,:,kFrame) , n1 , [] ) ...
                        ) , ...
                        outRes(1) , [] );
                end
            else
                for kFrame=1:nFrame
                    out(:,:,kFrame) = ...
                        reshape( ...
                        sum( ...
                        reshape( ...
                        reshape( ...
                        sum( ...
                        reshape( frame(:,:,kFrame) , n1 , [] ) ...
                        ) , ...
                        outRes(1) , [] ).' , ...
                        m2 , [] ) ...
                        ) , ...
                        outRes(2) , [] ).';
                end
            end
        end
        
        function out = polar3(theta,rho,z,varargin)
            %% POLAR3 Polar coordinate plot with color coded markers
            %
            % polar3(theta,rho,z) makes a plot using polar coordinates of
            % the angle THETA, in radians, versus the radius RHO. The color
            % of the markers is scaled according to the values in vector z.
            %
            % polar3(theta,rho,z,style) uses the marker specified in style
            %
            % polar3(...,'zMinMax',zBound) sets the z color scale limits to
            % the zBound values
            %
            % h = polar3(...) returns a handle to the plotted object in H.
            %
            % See also polar
            
            p = inputParser;
            p.addRequired('theta',@isnumeric);
            p.addRequired('rho',@isnumeric);
            p.addRequired('z',@isnumeric);
            p.addOptional('style','.',@ischar);
            p.addParamValue('zMinMax',[],@isnumeric);
            p.parse(theta,rho,z , varargin{:});
            style   = p.Results.style;
            zMinMax = p.Results.zMinMax;
            
            n = length(theta);
            if isempty(zMinMax)
                minZ = min(z);
                maxZ = max(z);
                fprintf(' @(utilities:polar3)> Z axis minmax: [%.2f,%.2f]\n',minZ,maxZ)
            else
                minZ = zMinMax(1);
                maxZ = zMinMax(2);
            end
            
            c    = colormap;
            nc = length(c);
            zc = fix((nc-1)*(z - minZ)/(maxZ-minZ) + 1);
            
            index = find(rho==max(rho));
            h = polar(theta(index),rho(index),'.');
            delete(h)
            
            h = zeros(n,1);
            hold on
            for k=1:n
                h(k) = polar(theta(k),rho(k),style);
                set(h(k),'zData',z(k),'color',c(zc(k),:))
            end
            hold off
            
            hc = colorbar;
            set(hc,'ylim',[minZ maxZ])
            set(get(hc,'children'),'YData',[minZ maxZ])
            
            if nargout == 1
                out = h;
            end
        end
        
        
        function out = defocusDistance(a4,focalLength,diameter,wavelength,unit)
            % DEFOCUSDISTANCE Focal point deplacement for a Zernike defocus
            %
            % out = defocusDistance(a4,focalLength,diameter,wavelength)
            % Compute the focal point relative position [meter] for the
            % Zernike (Noll normalized) focus coefficients [radian], the
            % focalLength [meter], the beam diameter [meter] and the
            % wavelength [meter] 
            %
            % out = defocusDistance(a4,focalLength,diameter,wavelength,unit) 
            % The result is converted into the appropriate unit: 3, 0, -3,
            % -6, -9 for example correspond to km, m, mm, micron, nm,
            % respectively
            
            out = 16*sqrt(3)*a4*(focalLength/diameter)^2/...
                ( 2*pi/wavelength - 16*sqrt(3)*a4*focalLength/diameter^2 );
            
            if nargin>4
                out = out*10^-unit;
            end
        end
        
        function out = outOfFocus(delta,focalLength,diameter,wavelength,unit)
            % OUTOFFOCUS Zernike focus for a focal point deplacement 
            %
            % out = outOfFocus(delta,focalLength,diameter,wavelength)
            % Compute the Zernike (Noll normalized) focus coefficients
            % [radian] for the focal point relative position [meter], the
            % focalLength [meter], the beam diameter [meter] and the
            % wavelength [meter]
            
            out = ( 2*pi*delta/wavelength ) / ...
                ( 16*sqrt(3)*( (focalLength/diameter)^2 + focalLength*delta/diameter^2 ) );
            
            if nargin>4
                out = (wavelength/(2*pi))*out*10^-unit;
            end
            
        end
        
        function out = orbitalVelocity(h,zen)
            %% ORBITALVELOCITY Orbital angular velocity
            %
            % out = orbitalVelocity(h) computes the orbital angular in
            % [rad/s] velocity at altitude h a zenith
            %
            % out = orbitalVelocity(h,zen) computes the orbital angular in
            % [rad/s] velocity at altitude h a zenith angle zen
            
            
            if nargin==1
                zen = 0;
            end
            out = sqrt(constants.G*constants.Me/(constants.Re+h)).*...
                (1-constants.Re*sin(zen)^2/(constants.Re+h))./h;
        end
        
        function out = pointAheadAngle(h,zen)
            %% POINTAHEADANGLE Point ahead angle
            %
            % out = pointAheadAngle(h) computes the orbital angular in
            % [rad] velocity at altitude h a zenith
            %
            % out = pointAheadAngle(h,zen) computes the orbital angular in
            % [rad] velocity at altitude h a zenith angle zen
            
            
            if nargin==1
                zen = 0;
            end
            out = 2*h*utilities.orbitalVelocity(h,zen)*sec(zen)/constants.c;
        end
        
        function [vertex,center,hp] = hexagonalArray(nCycle,pitch)
            %% HEXAGONALARRAY Array of hexagons
            %
            % [vertex,center] = hexagonalArray(nCycle,pitch) computes the
            % vertex and center coordinates of hexagons with the given
            % pitch arranged in a hexagonal array
            
            if nargin<2
                pitch=1;
            end
            a = pitch/sqrt(3);
            hexCoord = a*exp(1i*((0:5)*pi/3 + pi/2));
            count = 1;
            nSegment = 3*nCycle^2+3*nCycle+1;
            vertex = zeros(6,nSegment);
            vertex(:,count) = hexCoord;
            center = zeros(nSegment,1);
            for cycle=1:nCycle
                for o=1:6
                    zo = hexCoord + cycle*a*sqrt(3)*exp(1i*(o-1)*pi/3);
                    for k=1:cycle
                        zk = zo + (k-1)*a*sqrt(3)*exp(1i*((o-1)*pi/3+2*pi/3));
                        zk_center = mean(zk);
                        count = count + 1;
                        vertex(:,count) = zk;
                        center(count)   = zk_center;
                    end
                end
            end
            v = vertex(:);
            f = reshape(1:6*nSegment,6,nSegment);
            figure(nSegment)
            hp = patch('Faces',f','Vertices',[real(v(:)),imag(v(:))],'FaceColor',[1,1,1]*0.8);
%             line(real(center),imag(center),'color','r','marker','.')
            axis square
            set(gca,'ylim',get(gca,'xlim'))
            title(sprintf('%d segments',nSegment))
        end
        
        function B = eyeBlockDiag( A , n)
            %% EYEBLOCKDIAG Block diagonal concatenation
            %
            % B = eyeBlockDiag( A , n) concatenates n copies of the matrix A on
            % the diagonal of the matrix B. B is a sparse matrix.
            
            B = repmat( {sparse(A)} , 1 , n);
            B = blkdiag( B{:} ); 
        end
        
        function V = gramSchmidt(V)
            %% GRAMSCHMIDT Gram-Schmidt orthonormalization process
            %
            % V = gramSchmidt(V) orthonormalize the vector set V according
            % to the Gram-Schimdt process
            
            k = size(V,2);
            h = waitbar(0,'Gram-Schmith orthogonalization ...!');
            for j=1:k
                v = V(:,j);
                for i=1:j-1
                    u = V(:,i);
                    v = v - u*(u'*v)*(u'*u);
                end                
                V(:,j) = v/norm(v);
                waitbar(j/k)
            end
            close(h)
        end
        
        function out = besselJDerivative(nu,x)
            %% BESSELJDERIVATIVE Derivative of Bessel function of the first kind 
            %
            % out = besselJDerivative(nu,x) computes the derivative of the
            % Bessel function of the first kind of order n at x
            
            out = 0.5*( besselj(nu-1,x) - besselj(nu+1,x) );
            
        end
        
        function s = besselJDerivativeRoots(nu,ns)
           %% BESSELJDERIVATIVEROOUTS Roots of the drivative of Bessel function of the first kind 
           %
           % s = besselJDerivativeRoots(nu,ns) computes the ns first roots
           % of the derivative of the Bessel function of the first kind of
           % order n
            
            bjd = @(x) tools.besselJDerivative(nu,abs(x));
            s = zeros(1,ns);
            if nu==0
                s(1) = fzero( bjd , nu+3);
            else
                s(1) = fzero( bjd , nu);
            end
            for ks=2:ns
                x0 = ceil( s(ks-1) );
                bjd_x0 = bjd( x0 );
                x1 = x0 + pi; % the intervalle between 2 succesive roots tends towards pi but may larger for the first roots
                count = 0;
                while bjd_x0*bjd(x1)>0 && count<3
                    x1 = x1 + pi;
                    count = count + 1;
                end
                if count>2
                    error('Failed finding two points around roots where the signs of Bessel derivative differ!')
                end
                s(ks) = fzero( bjd , [x0,x1]);
            end
            
%             figure(103)
%             fplot( bjd, [0,ceil(s(end))])
%             grid
%             line( s, zeros(size(s)), 'color','r','marker','.','linestyle','none')

        end

        function rc = fitFwhm(profile)
            C = contourc(profile/max(profile(:)),[0.5,0.5]);
            rr = hypot(C(1,2:end),C(2,2:end));
            xc = sum(rr.*C(1,2:end))./sum(rr);
            yc = sum(rr.*C(2,2:end))./sum(rr);
%             line(xc,yc,'color','r','marker','x')
            rc = mean(sqrt((C(1,2:end)-xc).^2 + (C(2,2:end)-yc).^2));
        end
        
        function F = bilinearInterpolation(arg1,arg2,arg3,arg4,arg5)
            [nrows,ncols] = size(arg3);
            %     mx = numel(arg1); my = numel(arg2);
            s = 1 + (arg4-arg1(1))/(arg1(end)-arg1(1))*(ncols-1);
            t = 1 + (arg5-arg2(1))/(arg2(end)-arg2(1))*(nrows-1);
            
            
            % Matrix element indexing
            ndx = floor(t)+floor(s-1)*nrows;
            
            s(:) = (s - floor(s));
            t(:) = (t - floor(t));
            
            % Now interpolate.
            onemt = 1-t;
            F =  ( arg3(ndx).*(onemt) + arg3(ndx+1).*t ).*(1-s) + ...
                ( arg3(ndx+nrows).*(onemt) + arg3(ndx+(nrows+1)).*t ).*s;
        end
        
        function fr = gaussian(resolution,fwhm,n_f)
            
            u = (0:resolution-1)-(resolution)/2;
            [x,y] = meshgrid(u);
            r = hypot(x,y);
            sig = fwhm./(2.*sqrt(2*log(2)));
            f = exp(-r.^2/(2*sig.^2));
            f = f/sum(f(:));
            
%             n_f = 20;
            if nargin<3
                n_f = resolution;
            end
            if n_f<resolution/2
                fr = f;
                fr(:,end-n_f+1:end) =[];
                fr(end-n_f+1:end,:) =[];
                fr(:,1:n_f) =[];
                fr(1:n_f,:) =[];
                %     figure(101)
                %     subplot(1,2,1)
                %     imagesc(f)
                %     axis square
                %     subplot(1,2,2)
                %     imagesc(fr)
                %     axis square
            else
                fr = f;
            end
            
            
        end
        
        function out1 = pupAutoCorr(D,r)
            
            f_index       = r <= D;
            red         = r(f_index)./D;
            out1        = zeros(size(r));
            out1(f_index) = D.*D.*(acos(red)-red.*sqrt((1-red.*red)))./2;
            
        end
        
        function out2 = pupCrossCorr(R1,R2,r)
            
            out2 = zeros(size(r));
            f_index       = r <= abs(R1-R2);
            out2(f_index) = pi*min([R1,R2]).^2;
            
            f_index       = (r > abs(R1-R2)) & (r < (R1+R2));
            rho         = r(f_index);
            red         = (R1*R1-R2*R2+rho.*rho)./(2.*rho)/(R1);
            out2(f_index) = out2(f_index) + R1.*R1.*(acos(red)-red.*sqrt((1-red.*red)));
            red         = (R2*R2-R1*R1+rho.*rho)./(2.*rho)/(R2);
            out2(f_index) = out2(f_index) + R2.*R2.*(acos(red)-red.*sqrt((1-red.*red)));
            
        end
        
        function out = temporalSpectrum(nu,atm,spectrum)
            %% TEMPORALSPECTRUM Temporal power spectrum density
            %
            % out = phaseStats.temporalSpectrum(nu,spectrum) computes the
            % phase temporal power spectrum density from the spatial power
            % spectrum for a frozen flow atmosphere temporal model atm;
            % spectrum is a handle of an anonymous function fun(fr,fo,atm)
            
            out = zeros(size(nu));
            for kLayer = 1:atm.nLayer
                atmSlab = slab(atm,kLayer);
                [vx,vy] = pol2cart(atmSlab.layer.windDirection,atmSlab.layer.windSpeed);
                for k=1:numel(nu)
                    if vx>eps(atmSlab.layer.windSpeed)
                        out(k) = out(k) + quadgk( @integrandFy , -Inf, Inf);
                    else
                        out(k) = out(k) + quadgk( @integrandFx , -Inf, Inf);
                    end
                end
            end
            
            function int = integrandFy(fy)
                fx = (nu(k) -fy*vy)/vx;
                int = spectrum( hypot(fx,fy) , atan2(fy,fx), atmSlab)/vx;
            end
            
            function int = integrandFx(fx)
                fy = (nu(k) -fx*vx)/vy;
                int = spectrum( hypot(fx,fy) , atan2(fy,fx), atmSlab)/vy;
            end
        end
        
        function out = marechalStrehl(rmsWfe,band)
            
            out = (1 - (rmsWfe*2*pi/band.wavelength)^2/2)^2;
            
        end
        
        function [out,cvgce] = gerchbergSaxton(pupilPlaneImage,focalPlaneImage)
            
            source = sqrt(pupilPlaneImage);
            target = sqrt(focalPlaneImage);
            A = fftshift( ifft2( fftshift( target ) ) );
            phaseA = pi*(rand(size(source))*2-1);
%             figure,imagesc(abs(A))
            n = length(source);
            nIteration = 300;
            kIteration = 0;
            cvgce = zeros(1,nIteration);
            
            figure(111)
            subplot(2,4,[1,2])
            h(1) = imagesc(zeros(n,2*n));
            axis equal tight
            subplot(2,4,[5,6])
            h(2) = imagesc(zeros(n,2*n));
            axis equal tight
            subplot(2,4,[3,8])
            h(3) = imagesc(phaseA.*source);
            ht = title(sprintf('Iteration #: %d',0));
            axis equal tight
            colorbar('south')
            drawnow
            
            tic
            while kIteration<nIteration
                
                kIteration = kIteration + 1;
                
                B = source.*exp(1i*phaseA);
                C = fftshift( fft2( fftshift( B ) ) );
                D = target.*exp(1i*angle(C));
                A = fftshift( ifft2( fftshift( D ) ) );
                
                set(h(1),'CData',abs([C,D]))
                set(h(2),'CData',abs([B,A]))
                
                phaseA = angle(A);
                
                phaseA_ = phaseA;
                phaseA_(~source) = NaN;
                set(h(3),'CData',phaseA_)
                set(ht,'string',sprintf('Iteration #: %d',kIteration))
                drawnow
                
                cvgce(kIteration) = norm(abs(C).^2-focalPlaneImage,'fro');
                
            end
            toc
            
            out = phaseA;
        end
        
        function out = fftcentre(x,dir,n1,n2)
            %FFTCENTRE Computes the Fourier transform of a signal centered in the middle of the sample
            %The result is also centered on the same point. It is a conversion of an IDL routine written
            % by F. Cassaing (ONERA)
            %out = fftcentre(x,dir): x, the signal; dir: 1 for direct FT, -1 for inverse FT
            % IDL doc. written by Fred:
            % ;NOM :
            % ;   FFTCENTRE - Effectue une FFT d'un signal suppos? centr? sur le pixel m?dian
            % ;
            % ;CATEGORIE :
            % ;   Signal Processing Routines
            % ;
            % ;SYNTAXE :
            % ;   y=FFTCENTRE (x [,direction] [,/inverse] [,/double] [,/VERSION] [,/HELP])
            % ;
            % ;DESCRIPTION :
            % ;   Effectue une TF en  prenant pour origine non pas le pixel 0 mais le pixel
            % ;   central. Syntaxe identique ? FFT. Voir d?tails calcul ci-dessous.
            % ;
            % ;   ARGUMENTS :
            % ;          x : (entr?e) le signal ? transformer
            % ;  direction : (entr?e) le sens de la TF (-1 par d?faut=TF directe)
            % ;   /inverse : (entr?e) pour forcer direction ? +1
            % ;    /double : (entr?e) effectue le calcul en double
            % ; /overwrite : (entr?e) ?crase la variable x d'entr?e pour gagner en RAM
            % ;   /VERSION : (entr?e) affichage de la version avant l'ex?cution.
            % ;   /HELP    : (entr?e) affichage de la syntaxe et sortie du programme.
            % ;
            % ;   Principe:
            % ;   Soit une ?quence x(k) de longueur N points, centr?e sur #
            % ;      si N=2p    | 0 | 1 | 2 | . |p-1# p |p+1| . |2p-1|
            % ;      si N=2p+1  | 0 | 1 | 2 | . |p-1| # |p+1| . |2p-1|2p |
            % ;   Le pixel central tombe donc entre 2 pixels quand N est pair.
            % ;   Mais dans tous les cas, il a pour abscisse s=(N-1)/2.
            % ;
            % ;   La proc?dure est donc de recentrer le signal sur le pixel 0 avec un
            % ;   d?calage S1 vers la gauche de s=(N-1)/2 pixels, d'effectuer une FFT
            % ;   normale, puis enfin de recentrer la TF avec une translation S2 de (N-1)/2
            % ;   pixels vers la droite. La s?quence totale pour obtenir y=TF(x) est donc :
            % ;
            % ;        x(k)  --S1--> x'(k) --TF--> y'(l) --S2--> y(l)
            % ;
            % ;   Avec x(k+s)=x'(k) ou x'(k-s)=x(k) [idem pour y]. ON SUPPOSE DANS LA SUITE
            % ;   N IMPAIR DE SORTE QUE s SOIT ENTIER, MAIS CA MARCHE AUSSI POUR N PAIR.
            % ;
            % ;   La grosse ruse (merci Laurent Mugnier) et de remplacer le d?calage par une
            % ;   multiplication par un terme de phase dans l'espace conjug?. Ainsi, en
            % ;   appelant d (=-1 direct et =+1 inverse) la direction de la TF, et en notant
            % ;   w=exp[d2i pi/N] le terme de base du calcul de la TF, la s?quence
            % ;   pr?c?dente s'?crit :
            % ;
            % ;               y'(l)=SUM(k=0,N-1)   x'(k)  w^[  l   k]
            % ;     On translate k et l [muet] de s
            % ;             y'(l-s)=SUM(k=s,N-1+s) x'(k-s)w^[(l-s)(k-s)]
            % ;     On fait intervenir la def de x' et y'
            % ;                y(l)=SUM(k=s,N-1+s) x(k)   w^[(l-s)(k-s)]
            % ; Cette somme se coupe en 2 termes de k=+s ? N-1 et de k=N ? N-1+s. Par la
            % ; p?riodicit? des x(k) et des w^k, la derni?re se ram?ne ? k=0 ? s-1. D'o?
            % ;                y(l)=SUM(k=0,N-1)   x(k)   w^[(l-s)(k-s)]
            % ; ==> y(l)=p(-s)*p(l)*SUM(k=0,N-1) x(k)p(k) w^[  l    k]
            % ;
            % ;   o? p(k)=w^[-ks] est un phaseur. Il y a donc une rampe de phase p(k) ?
            % ;   appliquer avant la TF, la m?me rampe de phase p(l) ? appliquer apr?s la
            % ;   TF, ET un outsider, un terme de phase constant p(-s). Comme s=(N-1)/2, les
            % ;   deux rampes de phase s'?crivent exp[-dir*i*!pi*indice*(N-1)/N].
            % ;
            % ;   Le remplacement du d?calage par une multiplication permet peut-etre de
            % ;   gagner en temps d'ex?cution (? v?rifier) mais surtout de permettre le
            % ;   d?calage d'un nombre demi-entier de pixels. Je n'ai pas r?ussi ? d?montrer
            % ;   la formule pr?c?dente pour N pair (car avec s demi-entier et les deux
            % ;   s?quences x et x' ne d?coulent plus simplement l'une de l'autre) mais ?a
            % ;   marche num?riquement. Peut-?tre peut-on d?montrer en revenant ? des
            % ;   fonctions continues et en les num?risant ? ou faire jouer l'argument de la
            % ;   continuit?... Si quelqu'un trouve la d?monstration, elle sera la
            % ;   bienvenue!
            % ;
            % ;   Les 2 derni?res multiplications pourraient etre ?vit?e si l'on ne cherche
            % ;   que le module, ajouter un keyword ulterieur ?
            % ;
            % ;DIAGNOSTIC D'ERREUR :
            % ;
            % ;
            % ;VOIR AUSSI :
            % ;   FFT, la routine de base IDL
            % ;  VFFT, de L. Mugnier pour appeler FFTW (de L. Rousset-Rouvi?re) si besoin
            % ;
            % ;AUTEUR :
            % ;   $Author: cassaing $
            % ;
            % ;HISTORIQUE :
            % ;   $Log: fftcentre.pro,v $
            % ;   Revision 1.10  2001-02-02 16:43:01+01  cassaing
            % ;   Keyword double=double pass? a vfft supprim? car type variable pass? suffit
            % ;
            % ;   Revision 1.9  2001-01-26 14:55:22+01  cassaing
            % ;   Appel syst?matique de vfft au lieu de fft.
            % ;
            % ;   Revision 1.8  2001-01-17 12:29:54+01  cassaing
            % ;   Appel ? routine_courante() pour la doc
            % ;
            % ;   Revision 1.7  2001-01-17 11:18:35+01  cassaing
            % ;   Remplacement du message d'aide obsol?te par l'appel auto ? doc_library
            % ;
            % ;   Revision 1.6  2001-01-17 10:45:48+01  cassaing
            % ;   Phaseur 2D rempli par blas_axpy plutot que indgen#(fltarr+1). Bug doc corrig?
            % ;
            % ;   Revision 1.5  2001-01-16 17:59:00+01  cassaing
            % ;   VERSION PAS FINIE EN COURS DE DEBUG ....
            % ;
            % ;   Revision 1.4  2001-01-16 17:31:57+01  cassaing
            % ;   Prise en compte de /overwrite sur x (dej? fait ds fft sur var temp)
            % ;
            % ;   Revision 1.3  1999-09-20 10:35:51+02  cassaing
            % ;   Demo ?crite pour N impair, cas double/float group?s avec diripi
            % ;
            % ;   Revision 1.2  1999-09-17 16:57:28+02  cassaing
            % ;   Calcul enfin OK avec bons termes de phase. Mais d?mo th?orique pas claire...
            % ;
            % ;   Revision 1.1  1999-09-16 15:47:41+02  cassaing
            % ;   Initial revision
            % ;-
            
            %   R. Conan - LAOG-ONERA
            %   $Version 1.0 $  $Date: 2002/11/29 $
            
            persistent phasor coef
            
            if nargin<3
                [n1,n2,n3] = size(x);
            else
                x(n1,n2,size(x,3)) = 0;
                [n1,n2,n3] = size(x);
            end
            
            if isempty(phasor) || ndims(phasor)~=ndims(x) || any(size(phasor)~=size(x))
                
                disp('Info.: Pre-computing for fftcentre...')
                
                diripi = -dir.*1i.*pi;
                switch ndims(x)
                    case 1
                        n = length(x);
                        phasor = exp(diripi.*(0:(n-1)).*(1-n)./n);
                        coef   = exp(diripi*(n-1).^2./(2.*n));
                        %             switch dir
                        %                 case 1
                        %                     out    = coef.*fft(x.*phasor).*phasor;
                        %                 case -1
                        %                     out    = coef.*ifft(x.*phasor).*phasor;
                        %                 otherwise
                        %                     error('dir must be 1 or -1')
                        %             end
                    case 2
                        [u,v] = ndgrid((0:(n1-1)).*(1-n1)./n1,(0:(n2-1)).*(1-n2)./n2);
                        phasor = exp(diripi.*(u+v));
                        coef   = exp(diripi.*((n1-1).^2./n1+(n2-1).^2./n2)./2);
                        %             switch dir
                        %                 case 1
                        %                     out    = coef.*fft2(x.*phasor).*phasor;
                        %                 case -1
                        %                     out    = coef.*ifft2(x.*phasor).*phasor;
                        %                 otherwise
                        %                     error('dir must be 1 or -1')
                        %             end
                    case 3
                        [u,v] = ndgrid((0:(n1-1)).*(1-n1)./n1,(0:(n2-1)).*(1-n2)./n2);
                        phasor = exp(diripi.*(u+v));
                        coef   = exp(diripi.*((n1-1).^2./n1+(n2-1).^2./n2)./2);
                        phasor = repmat(phasor,[1,1,n3]);
                        %             switch dir
                        %                 case 1
                        %                     out    = coef.*fft2(x.*phasor).*phasor;
                        %                 case -1
                        %                     out    = coef.*ifft2(x.*phasor).*phasor;
                        %                 otherwise
                        %                     error('dir must be 1 or -1')
                        %             end
                    otherwise
                        error('fft dim must be < 3')
                end
                
            end
            
            switch ndims(x)
                case 1
                    switch dir
                        case 1
                            out    = coef.*fft(x.*phasor).*phasor;
                        case -1
                            out    = coef.*ifft(x.*phasor).*phasor;
                        otherwise
                            error('dir must be 1 or -1')
                    end
                case {2,3}
                    switch dir
                        case 1
                            out    = coef.*fft2(x.*phasor).*phasor;
                        case -1
                            out    = coef.*ifft2(x.*phasor).*phasor;
                        otherwise
                            error('dir must be 1 or -1')
                    end
                otherwise
                    error('fft dim must be < 3')
            end
                        
        end
        
        function out = oneParameterExample2(mu,alpha,q,p,a,nmax)
            
            n = 0:nmax;
            size_a = size(a);
            
            a = a(:)*0.5;
            
            coef1 = (-1).^n.*gamma( 2*(n+(alpha+mu)*0.5)/q ).*gamma( p-2*(n+(alpha+mu)*0.5)/q )./(gamma( n+1+alpha ).*gamma(n+1))./q;
            a1 = bsxfun( @power , a , 2*n+alpha+mu );
            a1 = bsxfun( @times, coef1, a1);
            
            coef2 = (-1).^n.*0.5.*gamma( -q*(n+p)*0.5+(alpha+mu)*0.5 ).*gamma( n+p )./(gamma( (alpha-mu)*0.5+1+q*(n+p)*0.5 ).*gamma(n+1));
            a2 = bsxfun( @power , a , q*(n+p) );
            a2 = bsxfun( @times, coef2, a2);
            
            out = sum( a1 + a2 , 2);
            out = 2.^mu.*out./gamma(p);
            out = reshape( out, size_a);
            
        end
        
        function out = oneParameterExample3(mu,alpha,q,p,a,nmax)
            
            n = 0:nmax;
            size_a = size(a);
            
            a = a(:);
            
            coef1 = (-1).^n.*gamma( 0.5+n+alpha ).*gamma( 2*(n+(2*alpha+mu)*0.5)/q ).*gamma( p-2*(n+(2*alpha+mu)*0.5)/q )./...
                (gamma( n+2*alpha+1 ).*gamma( n+alpha+1 ).*gamma(n+1))./q;
            a1 = bsxfun( @power , a , 2*n+2*alpha+mu );
            a1 = bsxfun( @times, coef1, a1);
            
            coef2 = (-1).^n.*0.5.*gamma( -q*(n+p)*0.5+(2*alpha+mu)*0.5 ).*gamma( 0.5+q*(n+p)*0.5-mu*0.5 ).*gamma( n+p )./...
                (gamma( (2*alpha-mu)*0.5+1+q*(n+p)*0.5 ).*gamma( 1+q*(n+p)*0.5-mu*0.5 ).*gamma(n+1));
            a2 = bsxfun( @power , a , q*(n+p) );
            a2 = bsxfun( @times, coef2, a2);
            
            out = sum( a1 + a2 , 2);
            out = out./(sqrt(pi)*gamma(p));
            out = reshape( out, size_a);
            
        end
        
    end

end
