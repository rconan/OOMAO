classdef gpuSource < source
    
    methods
        
        function obj = gpuSource(varargin)
            obj = obj@source(varargin{:});
            set(obj,'tag','GPU SOURCE');
        end
        
        function obj = mtimes(obj,otherObj)
            %% * Source object propagation operator
            %
            % src = src*otherObj propagate src through otherObj
            % multiplying the source amplitude by the otherObj transmitance
            % and adding the otherObj phase to the source phase
            
            nObj = numel(obj);
            for kObj = 1:nObj
                obj(kObj).opticalPath{ length(obj(kObj).opticalPath)+1 } = otherObj;
            end
            %             relay(otherObj,obj);
            fprintf(' *!* +source *!*\n')
        end
        
        function uplus(obj)
            %% UPLUS Source update/stream operator
            %
            % +obj streams the source object through its optical path
            
            nObj = numel(obj);
            for kObj = 1:nObj
                obj(kObj).mask        = [];
                obj(kObj).p_amplitude = 1;
                obj(kObj).p_phase     = 0;
            end
            %             cellfun(@(x)relay(x,obj),obj(1).opticalPath,'uniformOutput',false)
            gpuRun(obj)
        end
        
        function obj = times(obj,otherObj)
            %% .* Source object reset and propagation operator
            %
            % src = src.*otherObj propagates src through otherObj setting the
            % source amplitude to the otherObj transmitance and the source
            % phase to the otherObj phase
            
            mtimes(reset(obj),otherObj);
        end
        
    end
    
    methods (Access=private)
        
        function gpuRun(obj)
            
            nDevice = length( obj(1).opticalPath );
            nObj = size(obj,2);
            xL = zeros(nObj,1);
            yL = zeros(nObj,1);
            for kObj = 1:nObj
                xL(kObj) = obj(1,kObj,1).viewPoint(1);
                yL(kObj) = obj(1,kObj,1).viewPoint(2);
            end
            fr = cat(3,obj.extent);
            dmPhase = [];
            for kDevice = 1:nDevice
                deviceClass = class(obj(1).opticalPath{kDevice});
                switch deviceClass
                    case {'telescope','giantMagellanTelescope'}
                        tel = obj(1).opticalPath{kDevice};
%                         relayGpuSource(tel,obj)
                        set(obj,'mask',tel.pupilLogical)
                    case 'deformableMirror'
                        dm = obj(1).opticalPath{kDevice};
                        relay(dm,obj)
                        dmPhase = obj.phase;
                    case 'shackHartmann'
                        wfs = obj(1).opticalPath{kDevice};
                        nLenslet      = wfs.lenslets.nLenslet;
                        nValidLenslet = wfs.nValidLenslet;
                        nPixelLenslet = wfs.lenslets.nLensletImagePx;
                        nArray         = max( size(obj,2) , size(dmPhase,3) );
                        wfs.lenslets.nArray = nArray;
                        nLenslet2      = nLenslet^2;
                        nLenslets      = nLenslet^2*nArray;
                        nValidLenslets = nValidLenslet*nArray;
                        nLensletWavePx = nPixelLenslet;
                        nOutWavePx     = 2*nLensletWavePx + rem(nLensletWavePx,2);
                        n1             = nLensletWavePx*nLenslet;
                        n2             = n1*nArray;
                        
                        % phasor to align the spot on the center of the image
                        u = (0:(nLensletWavePx-1)).*(1-nLensletWavePx)./nOutWavePx;
                        phasor = exp(-1i*pi.*u);
                        phasor = phasor.'*phasor;
                        %%
                        index    = tools.rearrange( [n1,n1]  , [nLensletWavePx,nLensletWavePx] );
                        if isinf(tel.samplingTime)
                            error('oomao:gpuSource:gpuRun','Telescope sampling time must be finite!')
                        end
                        pupilWave = tel.pupil.*sqrt(tel.samplingTime*tel.area/tel.pixelArea)/nOutWavePx;
                        pupilWave = reshape( pupilWave(index) , nPixelLenslet,nPixelLenslet,nLenslet2);
                        if ~isempty(dmPhase)
                            fprintf('Reshaping the DM phase ...')
                            tId = tic;
                            buf2 = [];
                            for kDmPhase=1:size(dmPhase,3)
                                buf1 = dmPhase(:,:,kDmPhase);
                                buf1 = reshape( buf1(index) , nPixelLenslet,nPixelLenslet,nLenslet2);
                                buf1(:,:,~wfs.validLenslet) = [];
                                buf2 = cat(3,buf2,buf1);
                            end
                            dmPhase = gsingle(buf2);
                            clear('buf1','buf2')
                            et = toc(tId);
                            fprintf(' done in %4.2fs\n',et)
                        end
                        [xPup,yPup] = meshgrid(linspace(-1,1,tel.resolution)*tel.R);
                        yPup     = reshape( yPup(index)     , nPixelLenslet,nPixelLenslet,nLenslet2);
                        xPup     = reshape( xPup(index)     , nPixelLenslet,nPixelLenslet,nLenslet2);
                        buf2     = gzeros(nPixelLenslet,'single');
                        lensletIndex     = false(nOutWavePx);
                        lensletIndex(1:nPixelLenslet,1:nPixelLenslet) = true;
                        
                        %%
                        phasor           = gdouble( phasor);
                        pupilWave         = gsingle( pupilWave );
                        pupilWave = bsxfun( @times , pupilWave,phasor);
                        xPup             = gdouble( xPup );
                        yPup             = gdouble( yPup );
                        xL               = gdouble( xL );
                        yL               = gdouble( yL );
                        lensletIndex     = glogical( lensletIndex );
                        nHeight = size(obj,3);
                        if size(fr,3)==1
                            fr = repmat( fr , [1 ,1 ,nHeight] );
                        end
                        fr = gsingle(fr);
                        waveNumber = obj(1).waveNumber;
                        objectiveFocalLength = obj(1).objectiveFocalLength;
                        srcHeight  = gdouble( [obj(1,1,:).height]);
                        kLenslet_ = 1:nLenslets;
                        validLenslet = repmat( wfs.validLenslet(:)' , 1 , nArray );
                        kLenslet_( ~validLenslet ) = [];
                        kLenslet_ = gsingle( kLenslet_ );
                        naProfile = [obj(1,1,:).nPhoton];
                        [yLensletCoordinates,xLensletCoordinates] = ...
                            ndgrid(linspace(-1,1,nLenslet)*tel.R*(1-1/nLenslet));
                        xLensletCoordinates = gsingle(xLensletCoordinates);
                        yLensletCoordinates = gsingle(yLensletCoordinates);
                        %%
                        lensletIntensity = gzeros(nPixelLenslet*nPixelLenslet,nValidLenslets,'single');
                        nGpu = ceil(wfs.nValidLenslet*nOutWavePx^2*nArray*1e-6/135);
                        validLensletMaxStep = floor(nValidLenslets/nGpu);
                        validLensletRange   = 1:validLensletMaxStep;
                        
%                         F = gzeros(nPixelLenslet,nPixelLenslet,nValidLenslets);
                       
                        for kGpu=1:nGpu
                            
                            kValidLenslet_ = validLensletRange + (kGpu-1)*validLensletMaxStep;
                            kValidLenslet_(end) = min(kValidLenslet_(end),nValidLenslets);
                            %             gsync
                            fprintf(' -> GPU loop %d/%d: ',kGpu,nGpu)
                            tId = tic;
                            
                            if isempty(tel.opticalAberration) % WITHOUT ATMOSPHERE
                                
                                if isempty(fr)
                                    
                                    fprintf('(no transverse LGS profile) ')
                                    
                                    if isinf(obj.height) && isinf(obj.objectiveFocalLength) % NGS CASE
                                        
                                        fprintf('\b\b and NGS) ')
                                        
                                        gfor kValidLenslet = kValidLenslet_
                                        
%                                         for kHeight = 1:nHeight;
                                            
%                                             height          = srcHeight(kHeight);
%                                             naProfileHeight = naProfile(kHeight);
%                                             fr_             = fr(:,:,kHeight);
                                            
                                            kLenslet = kLenslet_(kValidLenslet);
                                            n = fix((kLenslet-1)/nLenslet2);
                                            k = kLenslet - n*nLenslet2;
                                            
%                                             lensletWave = pupilProp(xPup(:,:,k),yPup(:,:,k),xL(n+1),yL(n+1),...
%                                                 height,objectiveFocalLength,waveNumber,pupilWave(:,:,k));
                                            lensletWave = pupilWave(:,:,k).*exp(1i*dmPhase(:,:,kValidLenslet));
                                            
                                            buf1         = fft2( lensletWave , nOutWavePx , nOutWavePx );
%                                             buf2(:)      = abs( buf1(lensletIndex) ).^2;
                                            lensletIntensity(:,kValidLenslet) = abs( buf1(lensletIndex) ).^2;%lensletIntensity(:,kValidLenslet) + ...
                                                buf2(:);%*naProfileHeight;
                                            
%                                         end
                                        
                                        gend
                                        
                                    else
                                        
                                        gfor kValidLenslet = kValidLenslet_
                                        
                                        for kHeight = 1:nHeight;
                                            
                                            height          = srcHeight(kHeight);
                                            naProfileHeight = naProfile(kHeight);
                                            fr_             = fr(:,:,kHeight);
                                            
                                            kLenslet = kLenslet_(kValidLenslet);
                                            n = fix((kLenslet-1)/nLenslet2);
                                            k = kLenslet - n*nLenslet2;
                                            
                                            lensletWave = pupilProp(xPup(:,:,k),yPup(:,:,k),xL(n+1),yL(n+1),...
                                                height,objectiveFocalLength,waveNumber,pupilWave(:,:,k));
                                            
                                            buf1         = fft2( lensletWave , nOutWavePx , nOutWavePx );
                                            buf2(:)      = abs( buf1(lensletIndex) ).^2;
                                            lensletIntensity(:,kValidLenslet) = lensletIntensity(:,kValidLenslet) + ...
                                                buf2(:)*naProfileHeight;
                                            
                                        end
                                        
                                        gend
                                        
                                    end
                                    
                                else
                                    
                                    gfor kValidLenslet = kValidLenslet_
                                    
                                    for kHeight = 1:nHeight;
                                        
                                        height          = srcHeight(kHeight);
                                        naProfileHeight = naProfile(kHeight);
                                        fr_             = fr(:,:,kHeight);
                                        
                                        kLenslet = kLenslet_(kValidLenslet);
                                        n = fix((kLenslet-1)/nLenslet2);
                                        k = kLenslet - n*nLenslet2;
                                        
                                        lensletWave = pupilProp(xPup(:,:,k),yPup(:,:,k),xL(n+1),yL(n+1),...
                                            height,objectiveFocalLength,waveNumber,pupilWave(:,:,k));
                                        
                                        buf1         = fft2( lensletWave , nOutWavePx , nOutWavePx );
                                        buf2(:)      = abs( buf1(lensletIndex) ).^2;
                                        buf2         = conv2( buf2 , fr_ ,'same');
                                        lensletIntensity(:,kValidLenslet) = lensletIntensity(:,kValidLenslet) + ...
                                            buf2(:)*naProfileHeight;
                                        
                                    end
                                    
                                    gend
                                    
                                end
                                
                            else % WITH ATMOSPHERE
                                
                                fprintf('(with atmosphere ')
                                
                                atm_m           = tel.opticalAberration;
                                nLayer          = atm_m.nLayer;
                                layers          = atm_m.layer;
                                altitude_m      = gsingle([layers.altitude]);
                                sampler_m       = gsingle(linspace(-1,1,nPixelLenslet));
                                phase_m         = { layers.phase };
                                phase_m         = cellfun( @(x) gsingle(x) , phase_m , 'UniformOutput', false);
                                R_              = tel.R/nLenslet;
                                [ui,uj] = ndgrid(sampler_m*R_);
%                                 fprintf('R_=%4.2f) ',R_)
%                                 srcDirectionVector1 = src.directionVector(1);
%                                 srcDirectionVector2 = src.directionVector(2);
                                srcDirectionVector  = gsingle([obj(1,:,1).directionVector]);
%                                 srcHeight = src.height;
                                R_m    = gsingle([layers.D]*0.5);
                                nPixel = gsingle([layers.nPixel]);
                                
                                if isempty(fr)
                                    
                                    fprintf('but without transverse LGS profile) ')
                                    
                                    gfor kValidLenslet = kValidLenslet_
                                    
                                    for kHeight = 1:nHeight % src height
                                        
                                        height          = srcHeight(kHeight);
                                        naProfileHeight = naProfile(kHeight);
                                        fr_             = fr(:,:,kHeight);
                                        
                                        kLenslet = kLenslet_(kValidLenslet);
                                        n = fix((kLenslet-1)/nLenslet2);
                                        k = kLenslet - n*nLenslet2;
                                        
                                        xDvSrc = srcDirectionVector(1,n+1);
                                        yDvSrc = srcDirectionVector(2,n+1);
                                        yLC = yLensletCoordinates(k);
                                        xLC = xLensletCoordinates(k);
                                        
                                        for kLayer = 1:nLayer
                                            
                                            [ys,xs] = ndgrid(R_m(kLayer)*linspace(-1,1,nPixel(kLayer)));
                                            red = (1-altitude_m(kLayer)./height);
                                            xc = altitude_m(kLayer).*xDvSrc;
                                            yc = altitude_m(kLayer).*yDvSrc;
                                            xc = xc + yLC*red;
                                            yc = yc + xLC*red;
                                            yi = ui*red + xc;
                                            xi = uj*red + yc;
                                            
                                            arg3 = phase_m{kLayer};
                                            [nrows,ncols] = size(arg3);
                                            s = 1 + (xi-xs(1))/(xs(end)-xs(1))*(ncols-1);
                                            t = 1 + (yi-ys(1))/(ys(end)-ys(1))*(nrows-1);
                                            ndx = floor(t)+floor(s-1)*nrows;
                                            s(:) = (s - floor(s));
                                            t(:) = (t - floor(t));
                                            onemt = 1-t;
                                            F =  ( arg3(ndx).*(onemt) + arg3(ndx+1).*t ).*(1-s) + ...
                                                ( arg3(ndx+nrows).*(onemt) + arg3(ndx+(nrows+1)).*t ).*s;
                                            buf0 = pupilWave(:,:,k).*exp(1i*F);
                                            
                                            lensletWave = pupilProp(xPup(:,:,k),yPup(:,:,k),xL(n+1),yL(n+1),...
                                                height,objectiveFocalLength,waveNumber,buf0);
                                            
                                            buf1         = fft2( lensletWave , nOutWavePx , nOutWavePx );
                                            buf2(:)      = abs( buf1(lensletIndex) ).^2;
                                            lensletIntensity(:,kValidLenslet) = lensletIntensity(:,kValidLenslet) + ...
                                                buf2(:)*naProfileHeight;
                                            
                                        end
                                        
                                    end
                                    
                                    gend % gfor ends here!
                                    
                                else
                                    
                                    fprintf('and with transverse LGS profile) ')
                                    
                                    F = gzeros(nPixelLenslet,'single');

                                    gfor kValidLenslet = kValidLenslet_
                                    
                                    for kHeight = 1:nHeight % src height
                                        
                                        height          = srcHeight(kHeight);
                                        naProfileHeight = naProfile(kHeight);
                                        fr_             = fr(:,:,kHeight);
                                        
                                        kLenslet = kLenslet_(kValidLenslet);
                                        n = fix((kLenslet-1)/nLenslet2);
                                        k = kLenslet - n*nLenslet2;
                                        
                                        xDvSrc = srcDirectionVector(1,n+1);
                                        yDvSrc = srcDirectionVector(2,n+1);
                                        yLC = yLensletCoordinates(k);
                                        xLC = xLensletCoordinates(k);
                                        
                                        F(:) = 0;
                                        for kLayer = 1:nLayer
                                            
                                            [ys,xs] = ndgrid(R_m(kLayer)*linspace(-1,1,nPixel(kLayer)));
                                            red = (1-altitude_m(kLayer)./height);
                                            xc = altitude_m(kLayer).*xDvSrc;
                                            yc = altitude_m(kLayer).*yDvSrc;
                                            xc = xc + yLC*red;
                                            yc = yc + xLC*red;
                                            yi = ui*red + xc;
                                            xi = uj*red + yc;
                                            
                                            arg3 = phase_m{kLayer};
                                            [nrows,ncols] = size(arg3);
                                            s = 1 + (xi-xs(1))/(xs(end)-xs(1))*(ncols-1);
                                            t = 1 + (yi-ys(1))/(ys(end)-ys(1))*(nrows-1);
                                            ndx = floor(t)+floor(s-1)*nrows;
                                            s(:) = (s - floor(s));
                                            t(:) = (t - floor(t));
                                            onemt = 1-t;
                                            F =  F + ( arg3(ndx).*(onemt) + arg3(ndx+1).*t ).*(1-s) + ...
                                                ( arg3(ndx+nrows).*(onemt) + arg3(ndx+(nrows+1)).*t ).*s;
                                            
                                        end
                                        
                                        buf0 = pupilWave(:,:,k).*exp(1i*F);
                                        lensletWave = pupilProp(xPup(:,:,k),yPup(:,:,k),xL(n+1),yL(n+1),...
                                            height,objectiveFocalLength,waveNumber,buf0);
                                        
                                        buf1         = fft2( lensletWave , nOutWavePx , nOutWavePx );
                                        buf2(:)      = abs( buf1(lensletIndex) ).^2;
                                        buf2         = conv2( buf2 , fr_ ,'same');
                                        lensletIntensity(:,kValidLenslet) = lensletIntensity(:,kValidLenslet) + ...
                                            buf2(:)*naProfileHeight;
                                        
                                    end
                                    
                                    gend % gfor ends here!
                                    
                                end
                            end
                            
                            et = toc(tId);
                            fprintf('elapsed time: %4.2f\n',et)
                        end
%                         Fd = double(F);
%                         save('interpPhase.mat','Fd')
                        %                 geval(lensletIntensity)
                        %                 gsync
                        
                        %%
                        imagelets  = zeros(n1,n2);
                        index = tools.rearrange( [n1,n2] , [nLensletWavePx,nLensletWavePx] );
                        imagelets(index(:,validLenslet)) = wfs.lenslets.throughput*double(lensletIntensity(:));
                        wfs.lenslets.imagelets = imagelets;
                        
                        +wfs; %#ok<VUNUS>
                        
                end
            end
            
        end
        
    end
    
end

function lensletWave = pupilProp(xPup,yPup,xL,yL,...
    height,objectiveFocalLength,waveNumber,telPupil)

rho = hypot( xPup - xL , yPup - yL );
s  = hypot(rho,height);
s0 = hypot(rho,objectiveFocalLength);
sphericalPhase = waveNumber*(s - s0);

lensletWave  = telPupil.*exp(1i*sphericalPhase)./height;
end
