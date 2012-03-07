classdef gpuSource < source
    
    methods
        
        function obj = gpuSource(varargin)
            obj = obj@source(varargin{:});
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
            fr = obj(1).extent;
            for kDevice = 1:nDevice
                deviceClass = class(obj(1).opticalPath{kDevice});
                switch deviceClass
                    case {'telescope','giantMagellanTelescope'}
                        tel = obj(1).opticalPath{kDevice};
                    case 'shackHartmann'
                        wfs = obj(1).opticalPath{kDevice};
                end
            end
            nLenslet      = wfs.lenslets.nLenslet;
            nValidLenslet = wfs.nValidLenslet;
            nPixelLenslet = wfs.lenslets.nLensletImagePx;
            nArray         = size(obj,2);
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
            telPupil = tel.pupil.*sqrt(tel.samplingTime*tel.area/tel.pixelArea)/nOutWavePx;
            telPupil = reshape( telPupil(index) , nPixelLenslet,nPixelLenslet,nLenslet2);
            [xPup,yPup] = meshgrid(linspace(-1,1,tel.resolution)*tel.R);
            yPup     = reshape( yPup(index)     , nPixelLenslet,nPixelLenslet,nLenslet2);
            xPup     = reshape( xPup(index)     , nPixelLenslet,nPixelLenslet,nLenslet2);
            buf2     = gzeros(nPixelLenslet,'single');
            lensletIndex     = false(nOutWavePx);
            lensletIndex(1:nPixelLenslet,1:nPixelLenslet) = true;
            
            %%
            phasor           = gdouble( phasor);
            telPupil         = gsingle( telPupil );
            telPupil = bsxfun( @times , telPupil,phasor);
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
            srcHeight  = gdouble( [obj(1,1,:).height]/cos(tel.elevation) );
            kLenslet_ = 1:nLenslets;
            validLenslet = repmat( wfs.validLenslet(:)' , 1 , nArray );
            kLenslet_( ~validLenslet ) = [];
            kLenslet_ = gsingle( kLenslet_ );
            naProfile = [obj(1,1,:).nPhoton];
            %%
            lensletIntensity = gzeros(nPixelLenslet*nPixelLenslet,nValidLenslets,'single');
%             gsync
%             t = tic;
            for kHeight = 1:nHeight;
                
                height          = srcHeight(kHeight);
                naProfileHeight = naProfile(kHeight);
                fr_             = fr(:,:,kHeight);
                
                try
                    
                    gfor kValidLenslet = 1:nValidLenslets
                    
                    kLenslet = kLenslet_(kValidLenslet);
                    n = fix((kLenslet-1)/nLenslet2);
                    k = kLenslet - n*nLenslet2;
                    
                    lensletWave = pupilProp(xPup(:,:,k),yPup(:,:,k),xL(n+1),yL(n+1),...
                        height,objectiveFocalLength,waveNumber,telPupil(:,:,k));
                    
                    buf1         = fft2( lensletWave , nOutWavePx , nOutWavePx );
                    buf2(:)      = abs( buf1(lensletIndex) ).^2;
                    buf2         = conv2( buf2 , fr_ ,'same');
                    lensletIntensity(:,kValidLenslet) = lensletIntensity(:,kValidLenslet) + ...
                        buf2(:)*naProfileHeight;
                    
                    gend
                    %                 geval(lensletIntensity)
                    %                 gsync
                    
                catch err
                    
                    if ~strcmp(err.identifier,'jacket:runtime')
                        rethrow(err)
                    end
                    warning([err.identifier,': ',err.message])
                    
                    gfor kValidLenslet = 1:nValidLenslets/2
                    
                    kLenslet = kLenslet_(kValidLenslet);
                    n = fix((kLenslet-1)/nLenslet2);
                    k = kLenslet - n*nLenslet2;
                    
                    lensletWave = pupilProp(xPup(:,:,k),yPup(:,:,k),xL(n+1),yL(n+1),...
                        height,objectiveFocalLength,waveNumber,telPupil(:,:,k));
                    
                    buf1         = fft2( lensletWave , nOutWavePx , nOutWavePx );
                    buf2(:)      = abs( buf1(lensletIndex) ).^2;
                    buf2         = conv2( buf2 , fr_ ,'same');
                    lensletIntensity(:,kValidLenslet) = lensletIntensity(:,kValidLenslet) + ...
                        buf2(:)*naProfileHeight;
                    
                    gend
                    gfor kValidLenslet = 1+nValidLenslets/2:nValidLenslets
                    
                    kLenslet = kLenslet_(kValidLenslet);
                    n = fix((kLenslet-1)/nLenslet2);
                    k = kLenslet - n*nLenslet2;
                    
                    lensletWave = pupilProp(xPup(:,:,k),yPup(:,:,k),xL(n+1),yL(n+1),...
                        height,objectiveFocalLength,waveNumber,telPupil(:,:,k));
                    
                    buf1         = fft2( lensletWave , nOutWavePx , nOutWavePx );
                    buf2(:)      = abs( buf1(lensletIndex) ).^2;
                    buf2         = conv2( buf2 , fr_ ,'same');
                    lensletIntensity(:,kValidLenslet) = lensletIntensity(:,kValidLenslet) + ...
                        buf2(:)*naProfileHeight;
                    
                    gend
                    
                end
                
            end
%             toc(t)
            %%
            imagelets  = zeros(n1,n2);
            index = tools.rearrange( [n1,n2] , [nLensletWavePx,nLensletWavePx] );
            imagelets(index(:,validLenslet)) = wfs.lenslets.throughput*double(lensletIntensity(:));
            wfs.lenslets.imagelets = imagelets;
            
            +wfs; %#ok<VUNUS>
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
