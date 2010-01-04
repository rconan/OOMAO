%% shackHartmann Class Definition
% Shack-Hartmann wavefront sensor class
classdef lensletProcessing < handle

    properties
        % logbook id
        logsId=-1;
        % valid lenslet mask       
        validLenslet;
        % camera flat field
        flatField = 0;
        % camera pixel gains
        pixelGains = 1;
        % use quad-cell
        quadCell = false;
%         %% 
%         % * phase gradients
%         quadCellSlopes;
        % use center of gravity
        centroiding = true;
%         %% 
%         % * phase gradients
%         centroidingSlopes;
        % use matched filter
        matchedFilter = false;
%         %% 
%         % * phase gradients
%         matchedFilterSlopes;
        % use correlation
        correlation = false;
%         %% 
%         % * phase gradients
%         correlationSlopes;
%         %%
%         % * input wave
%         wave;
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
        % measurements reference
        referenceSlopes = 0;
    end
    
    properties (SetObservable=true)
        % measurements 
        slopes = 0;
    end
    
    properties (Dependent)
        % number of valid lenslet
        nValidLenslet;
        % intensity in each lenslet
        lensletIntensity;
    end
    
    properties (GetAccess=private,SetAccess=private)
        % index array to reshape a detector frame into a matrix with one
        % raster imagelet per column
        indexRasterLenslet = NaN;
        % lenslet centers
        lensletCenterX;
        lensletCenterY;
    end
    
    events
        logs;
    end

    methods
        
        % Constructor
        function obj = lensletProcessing(lenslets)
            error(nargchk(1, 3, nargin))
            % Centering the slopes at the lenslets center
            obj.referenceSlopes = ...
                ones(obj.nValidLenslet*2,1).*...
                (lenslets.nLensletWavePx-1)/2;
            % Slopes listener
            obj.slopesListener = addlistener(obj,'slopes','PostSet',...
                @(src,evnt) slopesDisplay(obj) );
            obj.slopesListener.Enabled = false;            
            % intensity listener
            obj.intensityListener = addlistener(obj.camera,'frame','PostSet',...
                @(src,evnt) intensityDisplay(obj) );
            obj.intensityListener.Enabled = false;            
            % Timer settings
            obj.paceMaker = timer;
            obj.paceMaker.name = 'Shack-Hartmann Wavefront Sensor';
            obj.paceMaker.TimerFcn = {@timerCallBack, obj};
            obj.paceMaker.ExecutionMode = 'FixedSpacing';
            obj.paceMaker.BusyMode = 'drop';
            obj.paceMaker.Period = 1e-1;
            obj.paceMaker.ErrorFcn = 'disp('' @detector: frame rate too high!'')';
            function timerCallBack( timerObj, event, a)
%                 fprintf(' @detector: %3.2fs\n',timerObj.instantPeriod)
                dataProcessing(a)
            end
            logBook.signIn(obj);
        end
       
        % Destructor
        function delete(obj)
            if isvalid(obj.paceMaker)
                if strcmp(obj.paceMaker.Running,'on')
                    stop(obj.paceMaker)
                end
                delete(obj.paceMaker)
            end
            if ishandle(obj.slopesDisplayHandle)
                delete(obj.slopesDisplayHandle)
            end
            if ishandle(obj.intensityDisplayHandle)
                delete(obj.intensityDisplayHandle)
            end
            logBook.add(obj,'Terminated!')
        end
        
        % Number of valid lenslet
        function nValidLenslet = get.nValidLenslet(obj)
            nValidLenslet = sum(obj.validLenslet(:));
        end
        
        % there is one
        function set.referenceSlopes(obj,val)
            obj.referenceSlopes = val;
            if ishandle(obj.slopesDisplayHandle)
                hc = get(obj.slopesDisplayHandle,'children');
                u = obj.referenceSlopes(1:end/2)+obj.lensletCenterX;
                v = obj.referenceSlopes(1+end/2:end)+obj.lensletCenterY;
                set(hc(2),'xData',u,'yData',v)
            end
        end
        
        % Computes the intensity in each lenslet
        function lensletIntensity = get.lensletIntensity(obj)
            if isempty(obj.camera.frame)
                lensletIntensity = [];
            else
                [nPx,mPx]  = size(obj.camera.frame);
                nPxLenslet = nPx/obj.lenslets.nLenslet;
                mPxLenslet = mPx/obj.lenslets.nLenslet;
                try
                    buffer     = obj.camera.frame(obj.indexRasterLenslet);
                catch ME
                    obj.indexRasterLenslet ...
                        = rearrange([nPx,mPx],[nPxLenslet,mPxLenslet]);
                    obj.indexRasterLenslet(:,~obj.validLenslet(:)) ...
                        = [];
                    buffer     = obj.camera.frame(obj.indexRasterLenslet);
                end
                lensletIntensity = sum(buffer);
            end
        end
        
        % Slopes computing
        function varargout = dataProcessing(obj)
            [nPx,mPx]  = size(obj.camera.frame);
            nPxLenslet = nPx/obj.lenslets.nLenslet;
            mPxLenslet = mPx/obj.lenslets.nLenslet;
            try
                buffer     = obj.camera.frame(obj.indexRasterLenslet);
            catch ME
                obj.indexRasterLenslet ...
                    = rearrange([nPx,mPx],[nPxLenslet,mPxLenslet]);
                obj.indexRasterLenslet(:,~obj.validLenslet(:)) ...
                    = [];
                buffer     = obj.camera.frame(obj.indexRasterLenslet);
            end
            % Buffer pre-processing
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
                index               = massLenslet~=0;
                massLenslet(~index) = [];
                buffer(:,~index)    = [];
                [x,y]               = ...
                    ndgrid((0:(nPxLenslet-1)),(0:(mPxLenslet-1)));
                xyBuffer  ...
                                    = zeros(2*obj.nValidLenslet,1);
                xBuffer             = bsxfun( @times , buffer , x(:) )  ;
                xBuffer             = sum( xBuffer ) ./ massLenslet  ;
                yBuffer             = bsxfun( @times , buffer , y(:) )  ;
                yBuffer             = sum( yBuffer ) ./ massLenslet  ;
                xyBuffer([index index]') ...
                                    = [xBuffer yBuffer]';
                obj.slopes          = xyBuffer - obj.referenceSlopes;     
            elseif obj.matchedFilter
            elseif obj.correlation
            end
            obj.slopes = obj.slopes.*obj.slopesUnits;
            if nargout>0
                varargout{1} = obj.slopes;
            end
        end

        %%
        % * propagate through the lenslet array, grab the detector frame
        % and process the frame
        function varargout = grabAndProcess(obj)
%             propagateThrough(obj.lenslets);
            grab(obj.camera)
            dataProcessing(obj);
            if nargout>0
                varargout{1} = obj.slopes;
            end
        end

        %% quiver like display of slopes
        function varargout = slopesDisplay(obj,varargin)
            if ishandle(obj.slopesDisplayHandle)
                if nargin>1
                    set(obj.slopesDisplayHandle,varargin{:})
                end
                hc = get(obj.slopesDisplayHandle,'children');
                u = obj.referenceSlopes(1:end/2)+...
                    obj.slopes(1:end/2)+obj.lensletCenterX;
                v = obj.referenceSlopes(1+end/2:end)+...
                    obj.slopes(1+end/2:end)+obj.lensletCenterY;
                set(hc(1),'xData',u,'yData',v)
            else
                obj.slopesDisplayHandle = hgtransform(varargin{:});
                [nPx,mPx]  = size(obj.camera.frame);
                nPxLenslet = nPx/obj.lenslets.nLenslet;
                mPxLenslet = mPx/obj.lenslets.nLenslet;
                % Display lenslet center with a cross
                u = (0:nPxLenslet:nPx-1);
                v = (0:mPxLenslet:mPx-1);
                [obj.lensletCenterX,obj.lensletCenterY] = ndgrid(u,v);
                obj.lensletCenterX = obj.lensletCenterX(obj.validLenslet(:));
                obj.lensletCenterY = obj.lensletCenterY(obj.validLenslet(:));
                line(obj.lensletCenterX + (nPxLenslet-1)/2,...
                    obj.lensletCenterY + (mPxLenslet-1)/2,...
                    'color','k','Marker','o',...
                    'linestyle','none',...
                    'parent',obj.slopesDisplayHandle)
                set(gca,'xlim',[0,nPx],...
                    'ylim',[0,mPx],'visible','off')
                axis square
                % Display slopes reference
                u = obj.referenceSlopes(1:end/2)+obj.lensletCenterX;
                v = obj.referenceSlopes(1+end/2:end)+obj.lensletCenterY;
                line(u,v,'color','b','marker','x',...
                    'linestyle','none',...
                    'parent',obj.slopesDisplayHandle)
                % Display slopes 
                u = obj.referenceSlopes(1:end/2)+...
                    obj.slopes(1:end/2)+obj.lensletCenterX;
                v = obj.referenceSlopes(1+end/2:end)+...
                    obj.slopes(1+end/2:end)+obj.lensletCenterY;
                line(u,v,'color','r','marker','+',...
                    'linestyle','none',...
                    'parent',obj.slopesDisplayHandle)
            end
            if nargout>0
                varargout{1} = obj.slopesDisplayHandle;
            end
        end
        
        function varargout = intensityDisplay(obj,varargin)
            intensity = zeros(obj.lenslets.nLenslet);
            intensity(obj.validLenslet) = obj.lensletIntensity;
            if ishandle(obj.intensityDisplayHandle)
                set(obj.intensityDisplayHandle,...
                    'Cdata',intensity,varargin{:})
            else
                obj.intensityDisplayHandle = imagesc(intensity,varargin{:});
                axis square xy
                colorbar
            end
            if nargout>0
                varargout{1} = obj.intensityDisplayHandle;
            end
        end
        
        function slopesAndFrameDisplay(obj,varargin)
            imagesc(obj.camera,varargin{:});
            slopesDisplay(obj,'matrix',makehgtform('translate',[1,1,0]),varargin{:});
        end
        
        function slopesAndIntensityDisplay(obj,varargin)
            intensityDisplay(obj,varargin{:});
            n  = obj.lenslets.nLensletImagePx;
            slopesDisplay(obj,'matrix',makehgtform('translate',-[(n-1)/2,(n-1)/2,0]/n,'scale',1/n,'translate',[1,1,0]*2),varargin{:});
        end
        
    end

end
