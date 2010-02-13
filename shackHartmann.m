classdef shackHartmann < handle
    % SHACKHARTMANN Create a Shack-Hartmann object
    %
    % obj = shackHartmann(nLenslet,detectorResolution) creates a
    % Shack-Hartmann object with a (nLenslet X nLenslet) lenslet array and
    % a detector with a detectorResolution resolution
    %
    % obj = shackHartmann(nLenslet,detectorResolution,validLenslet) creates
    % a Shack-Hartmann object with a (nLenslet X nLenslet) lenslet array, a
    % detector with a detectorResolution resolution and a logical mask of
    % size nLenslet setting the location of the valid lenslets inside the
    % lenslet array
    %
    % obj = shackHartmann(nLenslet,detectorResolution,validLenslet,guideStars) creates
    % a Shack-Hartmann object with a (nLenslet X nLenslet) lenslet array, a
    % detector with a detectorResolution resolution, a logical mask of
    % size nLenslet setting the location of the valid lenslets inside the
    % lenslet array and a guide star
    %
    % See also: lensletArray, detector, source, lensletArrayHowto,
    % detectorHowto
    
    properties
        % lenslet array object
        lenslets;
        % detector object
        camera;
        % camera flat field
        flatField = 0;
        % camera pixel gains
        pixelGains = 1;
        % use quad-cell
        quadCell = false;
        % use center of gravity
        centroiding = true;
        % use matched filter
        matchedFilter = false;
        % use correlation
        correlation = false;
        % timer
        paceMaker;
        % slopes display handle
        slopesDisplayHandle;
        % slope listener
        slopesListener;
        % intensity display handle
        intensityDisplayHandle;
        % intensity listener
        intensityListener;
        % frame pixel threshold
        framePixelThreshold = -inf;
        % slopes units (default:1 is pixel)
        slopesUnits = 1;
        % zernike to slopes conversion matrix
        zern2slopes;
        % wavefront sensor tag
        tag = 'SHACK-HARTMANN WAVEFRONT SENSOR';
    end
    
    properties (Dependent,SetObservable=true,GetObservable=true)
        % measurements
        slopes = 0;
    end
    
    properties (Dependent)
        % valid lenslet mask
        validLenslet;
        % measurements reference
        referenceSlopes;
    end
    
    properties (Dependent, SetAccess=private)
        % number of valid lenslet
        nValidLenslet;
        % number of slopes
        nSlope;
        % intensity in each lenslet
        lensletIntensity;
        % valid actuatord
        validActuator;
        % zernike coefficients
        zernCoefs;
    end
    
    properties (Access=private)
        p_slopes;
        p_referenceSlopes = 0;
        p_validLenslet;
        % index array to reshape a detector frame into a matrix with one
        % raster imagelet per column
        indexRasterLenslet = NaN;
        % lenslet centers
        lensletCenterX;
        lensletCenterY;
        log;
    end
    
    methods
        
        %% Constructor
        function obj = shackHartmann(nLenslet,detectorResolution,minLightRatio)
            error(nargchk(1, 4, nargin))
            obj.lenslets = lensletArray(nLenslet);
            obj.camera   = detector(detectorResolution);
            obj.lenslets.nLensletWavePx = ...
                detectorResolution/nLenslet;
            if nargin>2
                obj.lenslets.minLightRatio = minLightRatio;
            else
                obj.validLenslet = true(nLenslet);
            end
            obj.camera.frameGrabber ...
                = obj.lenslets;
            obj.referenceSlopes = zeros(obj.nValidLenslet*2,1);
            obj.p_referenceSlopes = ...
                repmat(obj.p_referenceSlopes,obj.lenslets.nArray,1);
            % Slopes listener
            obj.slopesListener = addlistener(obj,'slopes','PostSet',...
                @(src,evnt) slopesDisplay(obj) );
            obj.slopesListener.Enabled = false;
%             % intensity listener (BROKEN: shackhartmann is not deleted after a clear)
%             obj.intensityListener = addlistener(obj.camera,'frame','PostSet',...
%                 @(src,evnt) intensityDisplay(obj) );
% %             obj.intensityListener.Enabled = false;
            % Timer settings
            obj.paceMaker = timer;
            obj.paceMaker.name = 'Shack-Hartmann Wavefront Sensor';
%             obj.paceMaker.TimerFcn = {@timerCallBack, obj};(BROKEN: shackhartmann is not deleted after a clear)
            obj.paceMaker.ExecutionMode = 'FixedSpacing';
            obj.paceMaker.BusyMode = 'drop';
            obj.paceMaker.Period = 1e-1;
            obj.paceMaker.ErrorFcn = 'disp('' @detector: frame rate too high!'')';
            obj.log = logBook.checkIn(obj);
%             function timerCallBack( timerObj, event, a)
%                 %                 fprintf(' @detector: %3.2fs\n',timerObj.instantPeriod)
%                 a.grabAndProcess
%             end
        end
        
        %% Destructor
        function delete(obj)
%             if isvalid(obj.slopesListener)
%                 delete(obj.slopesListener)
%             end
%             if isvalid(obj.intensityListener)
%                 delete(obj.intensityListener)
%             end
%             if isvalid(obj.paceMaker)
%                 if strcmp(obj.paceMaker.Running,'on')
%                     stop(obj.paceMaker)
%                 end
%                 delete(obj.paceMaker)
%             end
            if ishandle(obj.slopesDisplayHandle)
                delete(obj.slopesDisplayHandle)
            end
            if ishandle(obj.intensityDisplayHandle)
                delete(obj.intensityDisplayHandle)
            end
            delete(obj.lenslets)
            delete(obj.camera)
            checkOut(obj.log,obj)
        end
        
        function display(obj)
            %% DISPLAY Display object information
            %
            % disp(obj) prints information about the Shack-Hartmann
            % wavefront sensor object
          
            fprintf('___ %s ___\n',obj.tag)
            fprintf(' Shack-Hartmann wavefront sensor: \n  . %d lenslets total on the pupil\n  . %d pixels per lenslet \n',...
                obj.nValidLenslet,obj.camera.resolution(1)/obj.lenslets.nLenslet)
            if isinf(obj.framePixelThreshold)
                algoProp = ', no thresholding!';
            else
                algoProp = sprintf(', pixel threshold: %d\n',obj.framePixelThreshold);
            end
            
            algo = {'quadCell','centroiding','matchedFilter','correlation'};
            algoTF = [obj.quadCell,obj.centroiding,obj.matchedFilter,obj.correlation];
            fprintf('  . spot algorithm: %s%s\n',algo{algoTF},algoProp);
            fprintf('----------------------------------------------------\n')
            display(obj.lenslets)
            display(obj.camera)

        end
        
        %% Get and Set slopes
        function slopes = get.slopes(obj)
            slopes = obj.p_slopes;
        end
        function set.slopes(obj,val)
            obj.p_slopes = val;
        end
        
        %% Get and Set valid lenslets
        function validLenslet = get.validLenslet(obj)
            validLenslet = obj.p_validLenslet;
        end
        function set.validLenslet(obj,val)
            obj.p_validLenslet = val;
            obj.referenceSlopes(~[obj.validLenslet(:);obj.validLenslet(:)]) = [];
        end
        
        %% Get number of valid lenslet
        function nValidLenslet = get.nValidLenslet(obj)
            nValidLenslet = sum(obj.validLenslet(:));
        end
        
        %% Get number of slopes
        function nSlope = get.nSlope(obj)
            nSlope = obj.nValidLenslet*2;
        end
        
        %% Get valid actuators
        function val = get.validActuator(obj)
            nElements            = 2*obj.lenslets.nLenslet+1; % Linear number of lenslet+actuator
            validLensletActuator = zeros(nElements);
            index                = 2:2:nElements; % Lenslet index
            validLensletActuator(index,index) = obj.validLenslet;
            for xLenslet = index
                for yLenslet = index
                    if validLensletActuator(xLenslet,yLenslet)==1
                        xActuatorIndice = [xLenslet-1,xLenslet-1,...
                            xLenslet+1,xLenslet+1];
                        yActuatorIndice = [yLenslet-1,yLenslet+1,...
                            yLenslet+1,yLenslet-1];
                        validLensletActuator(xActuatorIndice,yActuatorIndice) = 1;
                    end
                end
            end
            index = 1:2:nElements; % Actuator index
            val   = logical(validLensletActuator(index,index));
        end
        
        %% Get/Set the reference spots and update spots location display if
        %% there is one
        function val = get.referenceSlopes(obj)
            val = obj.p_referenceSlopes;
        end
        function set.referenceSlopes(obj,val)
            obj.p_referenceSlopes = val;
            if ishandle(obj.slopesDisplayHandle)
                hc = get(obj.slopesDisplayHandle,'children');
                u = obj.p_referenceSlopes(1:end/2)+obj.lensletCenterX;
                v = obj.p_referenceSlopes(1+end/2:end)+obj.lensletCenterY;
                set(hc(2),'xData',u,'yData',v)
            end
        end
        
        % Get the zernike coeficients
        function val = get.zernCoefs(obj)
            val = obj.zern2slopes\obj.slopes;
            val(1,:) = []; % piston=0 removed
        end
        
        %% Computes the intensity in each lenslet
        function lensletIntensity = get.lensletIntensity(obj)
            if isempty(obj.camera.frame)
                lensletIntensity = [];
            else
                [nPx,mPx]  = size(obj.camera.frame);
                nLensletArray = obj.lenslets.nArray;
                nPxLenslet = nPx/obj.lenslets.nLenslet;
                mPxLenslet = mPx/obj.lenslets.nLenslet/nLensletArray;
                try
                    buffer     = obj.camera.frame(obj.indexRasterLenslet);
                catch ME
                    fprintf( '@(shackHartmann)> %s\n',ME.identifier)
                    obj.indexRasterLenslet ...
                        = utilities.rearrange([nPx,mPx],[nPxLenslet,mPxLenslet]);
                    v = ~obj.validLenslet(:);
                    v = repmat(v,nLensletArray,1);
                    obj.indexRasterLenslet(:,v) = [];
                    buffer     = obj.camera.frame(obj.indexRasterLenslet);
                end
                lensletIntensity = sum(buffer);
            end
        end
        
        
        function setValidLenslet(obj,pupilIntensity)
            %% SETVALIDLENSLET Valid lenslet mask
            %
            % setValidLenslet(obj,pupilIntensity) sets the mask of valid
            % lenslet based on the value of minLightRatio in the lenslets
            % object providing the pupil intensity map
            
            n = length(pupilIntensity);
            nL = n/obj.lenslets.nLenslet;
            pupilIntensity = reshape(pupilIntensity,nL,n*obj.lenslets.nLenslet);
            pupilIntensity = sum(pupilIntensity);
            pupilIntensity = reshape(pupilIntensity,obj.lenslets.nLenslet,obj.lenslets.nLenslet*nL);
            pupilIntensity = reshape(pupilIntensity',nL,obj.lenslets.nLenslet^2);
            pupilIntensity = sum(pupilIntensity);
            surfPx         = nL^2;
            pupilIntensity = pupilIntensity/surfPx;
            obj.p_validLenslet  = logical( ...
                reshape( pupilIntensity>=obj.lenslets.minLightRatio , ...
                obj.lenslets.nLenslet,obj.lenslets.nLenslet));
            obj.referenceSlopes = zeros(2*obj.nValidLenslet,1);
            obj.p_referenceSlopes = ...
                repmat(obj.p_referenceSlopes,obj.lenslets.nArray,1);
        end

        
        function varargout = dataProcessing(obj)
            %% DATAPROCESSING Processing a SH-WFS detector frame
            %
            % dataProcessing(obj) computes the WFS slopes
            %
            % out = dataProcessing(obj) computes and returns the WFS slopes
            
            [nPx,mPx,nFrame]  = size(obj.camera.frame);
            nLensletArray = obj.lenslets.nArray;
            nPxLenslet = nPx/obj.lenslets.nLenslet;
            mPxLenslet = mPx/obj.lenslets.nLenslet/nLensletArray;
            if numel(obj.indexRasterLenslet)~=(nPxLenslet*mPxLenslet*obj.nValidLenslet*nLensletArray*nFrame)
                %             try
                % %                 u = obj.indexRasterLenslet;
                % %                 if nFrame>1
                % %                     u = repmat(u,[1,1,nFrame]);
                % %                 end
                %                 buffer     = obj.camera.frame(obj.indexRasterLenslet);
                %             catch ME
                fprintf( '@(shackHartmann)> Setting the raster index \n')
                % get lenslet index
                obj.indexRasterLenslet ...
                    = utilities.rearrange([nPx,mPx/nLensletArray,nLensletArray*nFrame],[nPxLenslet,mPxLenslet]);
                % remove index from non-valid lenslets
                v = ~obj.validLenslet(:);
                v = repmat(v,nLensletArray,1);
                v = repmat(v,nFrame,1);
                obj.indexRasterLenslet(:,v) = [];
                %                 u = obj.indexRasterLenslet;
                %                 if nFrame>1
                %                     u = repmat(u,[1,1,nFrame]);
                %                 end
            end
            % Buffer pre-processing
            buffer     = obj.camera.frame(obj.indexRasterLenslet);
            buffer = (buffer - obj.flatField)./obj.pixelGains;
            % Thresholding
            if isfinite(obj.framePixelThreshold)
                buffer           = buffer - obj.framePixelThreshold;
                buffer(buffer<0) = 0;
            end
            % Centroiding
            if obj.quadCell
            elseif obj.centroiding
                massLenslet         = sum(buffer);
                %                 massLenslet(~index) = [];
                %                 buffer(:,~index)    = [];
                %                 size(buffer)
                [x,y]               = ...
                    meshgrid((0:(nPxLenslet-1)),(0:(mPxLenslet-1)));
                %                 xyBuffer  ...
                %                     = zeros(2*obj.nValidLenslet,1);
                xBuffer             = bsxfun( @times , buffer , x(:) )  ;
                xBuffer             = sum( xBuffer ) ./ massLenslet  ;
%                 xBuffer             = squeeze((xBuffer));
                yBuffer             = bsxfun( @times , buffer , y(:) )  ;
                yBuffer             = sum( yBuffer ) ./ massLenslet  ;
%                 yBuffer             = squeeze((yBuffer));
%                 xyBuffer = squeeze([xBuffer  yBuffer]);
%                 size(xyBuffer)
                xBuffer = reshape(xBuffer,obj.nValidLenslet,nLensletArray*nFrame);
                yBuffer = reshape(yBuffer,obj.nValidLenslet,nLensletArray*nFrame);
                sBuffer = bsxfun(@minus,[xBuffer ; yBuffer],obj.referenceSlopes).*obj.slopesUnits;
                index = isnan(sBuffer);
                if any(index(:)) % if all pixels threshold
                    warning('OOMAO:shackHartmann:dataProcessing',...
                        'Threshold is probably too high1')
                    sBuffer(index) = obj.p_slopes(index);
                end
                obj.p_slopes = sBuffer;
            elseif obj.matchedFilter
            elseif obj.correlation
            end
            if nargout>0
                varargout{1} = obj.p_slopes;
            end
        end
        
        function zern = getZernike(obj,radialOrder)
            zern = zernike(1:zernike.nModeFromRadialOrder(radialOrder),...
                'resolution',obj.lenslets.nLenslet,...
                'pupil',double(obj.validLenslet));
            dzdxy = [zern.xDerivative(obj.validLenslet,:);zern.yDerivative(obj.validLenslet,:)];
            zern.c = dzdxy\obj.p_slopes;
%             zern.c = dzdxy'*obj.p_slopes;
        end
        
        function varargout = grabAndProcess(obj)
            %% GRABANDPROCESS Frame grabbing and processing
            %
            % grabAndProcess(obj) grabs a frame and computes the slopes
            %
            % out = grabAndProcess(obj) grabs a frame, computes and returns
            % the slopes
            
            warning('OOMAO;shackHartmann:grabAndProcess',...
                'DEPRECATED! Instead use the uplus operator (+obj)')
            grab(obj.camera)
            dataProcessing(obj);
            if nargout>0
                varargout{1} = obj.p_slopes;
            end
        end
        function varargout = uplus(obj)
            %% UPLUS + Update operator
            %
            % +obj grabs a frame and computes the slopes
            %
            % obj = +obj returns the shackHartmann object
            
            grab(obj.camera)
            dataProcessing(obj);
            if nargout>0
                varargout{1} = obj;
            end
        end
  
        function relay(obj,src)
            %% RELAY
            propagateThrough(obj.lenslets,src)
            grabAndProcess(obj)
        end
        
        function varargout = slopesDisplay(obj,varargin)
            %% SLOPESDISPLAY WFS slopes display
            %
            % slopesDisplay(obj) displays quiver plot of the slopes
            %
            % slopesDisplay(obj,'PropertyName',PropertyValue) displays
            % quiver plot of the slopes and set the properties of the
            % graphics object quiver
            %
            % h = slopesDisplay(obj,...) returns the graphics handle
            %
            % See also: quiver
            
            if ishandle(obj.slopesDisplayHandle)
                if nargin>1
                    set(obj.slopesDisplayHandle,varargin{:})
                end
                hc = get(obj.slopesDisplayHandle,'children');
                u = obj.referenceSlopes(1:end/2)+...
                    obj.p_slopes(1:end/2)+obj.lensletCenterX;
                v = obj.referenceSlopes(1+end/2:end)+...
                    obj.p_slopes(1+end/2:end)+obj.lensletCenterY;
                set(hc(1),'xData',u,'yData',v)
            else
                obj.slopesDisplayHandle = hgtransform(varargin{:});
                [nPx,mPx]  = size(obj.camera.frame);
                nLensletArray = obj.lenslets.nArray;
                nPxLenslet = nPx/obj.lenslets.nLenslet;
                mPxLenslet = mPx/obj.lenslets.nLenslet/nLensletArray;
                % Display lenslet center with a cross
                u = (0:nPxLenslet:nPx-1);
                v = (0:mPxLenslet:mPx-1);
                [obj.lensletCenterX,obj.lensletCenterY] = ndgrid(u,v);
                obj.lensletCenterX = obj.lensletCenterX(obj.validLenslet(:));
                obj.lensletCenterY = obj.lensletCenterY(obj.validLenslet(:));
                if nLensletArray>1
                    offset = (0:nLensletArray-1).*nPx;
                    obj.lensletCenterX = repmat(obj.lensletCenterX,1,nLensletArray);
                    obj.lensletCenterX = bsxfun(@plus,obj.lensletCenterX,offset);
                    obj.lensletCenterX = obj.lensletCenterX(:);
                    obj.lensletCenterY = repmat(obj.lensletCenterY,nLensletArray,1);
                end
                line(obj.lensletCenterX + (nPxLenslet-1)/2,...
                    obj.lensletCenterY + (mPxLenslet-1)/2,...
                    'color','k','Marker','.',...
                    'linestyle','none',...
                    'parent',obj.slopesDisplayHandle)
                set(gca,'xlim',[0,nPx],...
                    'ylim',[0,mPx*nLensletArray],'visible','off')
                axis equal tight
                % Display slopes reference
                u = obj.referenceSlopes(1:end/2)+obj.lensletCenterX;
                v = obj.referenceSlopes(1+end/2:end)+obj.lensletCenterY;
                line(u,v,'color','b','marker','x',...
                    'linestyle','none',...
                    'parent',obj.slopesDisplayHandle)
                % Display slopes
                u = obj.referenceSlopes(1:end/2)+...
                    obj.p_slopes(1:end/2)+obj.lensletCenterX;
                v = obj.referenceSlopes(1+end/2:end)+...
                    obj.p_slopes(1+end/2:end)+obj.lensletCenterY;
                line(u,v,'color','r','marker','+',...
                    'linestyle','none',...
                    'parent',obj.slopesDisplayHandle)
            end
            if nargout>0
                varargout{1} = obj.slopesDisplayHandle;
            end
        end
        
        function varargout = intensityDisplay(obj,varargin)
            %% INTENSITYDISPLAY WFS lenslet intensity display
            %
            % intensityDisplay(obj) displays the intensity of the lenslets
            %
            % intensityDisplay(obj,'PropertyName',PropertyValue) displays
            % the intensity of the lenslets and set the properties of the
            % graphics object imagesc
            %
            % h = intensityDisplay(obj,...) returns the graphics handle
            %
            % See also: imagesc
            
            intensity = nan(obj.lenslets.nLenslet,obj.lenslets.nLenslet*obj.lenslets.nArray);
            v = obj.validLenslet(:);
            v = repmat(v,obj.lenslets.nArray,1);
            intensity(v) = obj.lensletIntensity;
            if ishandle(obj.intensityDisplayHandle)
                set(obj.intensityDisplayHandle,...
                    'Cdata',intensity,varargin{:})
            else
                obj.intensityDisplayHandle = imagesc(intensity,varargin{:});
                axis equal tight xy
%                 set(gca,'Clim',[floor(min(intensity(v))),max(intensity(v))])
                colorbar
            end
            if nargout>0
                varargout{1} = obj.intensityDisplayHandle;
            end
        end
        
        function slopesAndFrameDisplay(obj,varargin)
            imagesc(obj.camera,varargin{:});
            slopesDisplay(obj,'matrix',...
                makehgtform('translate',[1,1,0]),varargin{:});
%             if isinf(obj.framePixelThreshold)
%                 clim = get(gca,'clim');
%                 set(gca,'clim',[obj.framePixelThreshold,clim(2)])
%             end
                
        end
        
        function slopesAndIntensityDisplay(obj,varargin)
            intensityDisplay(obj,varargin{:});
            n  = obj.lenslets.nLensletImagePx;
            slopesDisplay(obj,'matrix',...
                makehgtform('translate',-[(n-1)/2,(n-1)/2,0]/n,'scale',1/n,'translate',[1,1,0]*2),varargin{:});
        end
        
    end
    
end
