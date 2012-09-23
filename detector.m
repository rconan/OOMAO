classdef detector < handle
    % DETECTOR Create a detector object
    %
    % obj = detector(resolution) creates a detector object from the detector
    % resolution
    
    properties
        % Add or not photon noise to the frame        
        photonNoise = false;
        % Photon background noise (# of photon per frame)
        nPhotonBackground = 0;
        % # of photo-electron per pixel rms
        readOutNoise = 0;
        % quantum efficiency
        quantumEfficiency = 1;
        % units of one pixel
        pixelScale;
        % detector resolution
        resolution;
        % detector region of interest
        regionOfInterest;
        roiSouthWestCorner;
        % frame rate [Hz]
        frameRate = 1;
        % clock rate [Hz]
        clockRate = 1;
        %  Delay after which the camera start integrating
        startDelay = 0;
        % detector timer
        paceMaker;
        % frame update listener
        frameListener;
        % frame grabber callback function
        frameGrabber;
        % detector tag
        tag= 'DETECTOR';
        frameBuffer = 0;
        frameCount  = 0;
    end
    
    properties (Dependent)
        % Exposure time [second]
        exposureTime;
    end
    
    properties (SetObservable=true)
        % detector frame
        frame;
    end
    
    properties (Access=protected)
        frameHandle;
        log;
    end
    
    properties (Access=private)
        p_exposureTime;
    end
    
    methods
        
        %% Constructor
        function obj = detector(resolution,pixelScale)
            error(nargchk(1, 2, nargin))
            if numel(resolution)==1
                obj.resolution = resolution*ones(1,2);
            else
                obj.resolution = resolution;
            end
            obj.regionOfInterest   = obj.resolution;
            obj.roiSouthWestCorner = [1,1];
            if nargin>1
                obj.pixelScale  = pixelScale;
            end
            obj.exposureTime = 1;

            setFrameListener(obj)
            
            % Timer settings
            obj.paceMaker = timer;
            obj.paceMaker.name = 'Detector';
%             obj.paceMaker.TimerFcn = @(src,evnt) obj.grab;% {@timerCallBack, obj};
            obj.paceMaker.ExecutionMode = 'FixedSpacing';
            %             obj.paceMaker.BusyMode = 'error';
            obj.paceMaker.Period = 1;
            obj.paceMaker.ErrorFcn = 'disp('' @detector: frame rate too high!'')';
%             function timerCallBack( timerObj, event, a)
%                 %                 fprintf(' @detector: %3.2fs\n',timerObj.instantPeriod)
%                 a.grab;
%             end
            %             obj.frameRate = 1;
            obj.log = logBook.checkIn(obj);
        end
        
        %% Destructor
        function delete(obj)
            if ~isempty(obj.frameHandle) && ishandle(obj.frameHandle(1))
                delete(get(obj.frameHandle(1),'Parent'));
            end
            if isvalid(obj.paceMaker)
                if strcmp(obj.paceMaker.Running,'on')
                    stop(obj.paceMaker)
                end
                delete(obj.paceMaker)
            end
            checkOut(obj.log,obj)
        end
        
        function display(obj)
            %% DISPLAY Display object information
            %
            % disp(obj) prints information about the detector object
          
            fprintf('___ %s ___\n',obj.tag)
            fprintf(' %dx%d pixels camera \n',...
                obj.resolution)
            if ~isempty(obj.pixelScale)
            fprintf('  . pixel scale: %4.2f milli-arcsec \n',...
                obj.pixelScale*constants.radian2arcsec*1000)                
            end            
            fprintf('  . quantum efficiency: %3.1f \n',...
                obj.quantumEfficiency)
            if obj.photonNoise
                fprintf('  . photon noise enabled\n')
            else
                fprintf('  . photon noise disabled\n')
            end
            fprintf('  . %.1f photo-events rms read-out noise \n',...
                obj.readOutNoise)
            fprintf('  . %3.1fms exposure time and %3.1fHz frame rate \n',...
                obj.exposureTime*1e3,obj.frameRate)
            fprintf('----------------------------------------------------\n')
            
        end
        
        function obj = saveobj(obj)
            %% SAVEOBJ
            delete(obj.frameListener)
            add(obj.log,obj,'Save!')
        end                
        
        %% Set/Get exposureTime
        function set.exposureTime(obj,val)
            obj.p_exposureTime = val;
            if obj.clockRate==1
                obj.clockRate = 1/obj.p_exposureTime;
                fprintf('Clock rate is: %.2fHz\n',obj.clockRate)
            end
        end
        function val = get.exposureTime(obj)
            val = obj.p_exposureTime;
        end
        
        function imagesc(obj,varargin)
            %% IMAGESC Display the detector frame
            %
            % imagesc(obj) displays the frame of the detector object
            %
            % imagesc(obj,'PropertyName',PropertyValue) displays the frame of
            % the detector object and set the properties of the graphics object
            % imagesc
            %
            % h = imagesc(obj,...) returns the graphics handle
            %
            % See also: imagesc
            
            if ndims(obj.frame)<3 % frame must be a 2-D array for plotting
                
                if isempty(obj.frame)
                    m_frame = obj.frameBuffer;
                    titleColor = 'r';
                else
                    m_frame = obj.frame;
                    titleColor = 'k';
                end
                
                if all(ishandle(obj.frameHandle)) && ~isempty(obj.frameHandle)
                    set(obj.frameHandle(1),'Cdata',m_frame,varargin{:});
                    %                 xAxisLim = [0,size(obj.frame,2)]+0.5;
                    %                 yAxisLim = [0,size(obj.frame,1)]+0.5;
                    %                 set( get(obj.frameHandle,'parent') , ...
                    %                     'xlim',xAxisLim,'ylim',yAxisLim);
                    set(obj.frameHandle(2),'String',...
                        sprintf('Frame #%d',obj.frameCount),...
                        'color',titleColor)
                    [n,m] = size(m_frame);
                    set(get(obj.frameHandle(1),'parent'),'xlim',0.5+[0 m],'ylim',0.5+[0 n])
                else
                    if isempty(obj.pixelScale)
                        obj.frameHandle(1) = image(m_frame,...
                            'CDataMApping','Scaled',...
                            varargin{:});
                    else
                        [n,m] = size(m_frame);
                        x = 0.5*linspace(-1,1,n)*...
                            n*obj.pixelScale*constants.radian2arcsec;
                        y = 0.5*linspace(-1,1,m)*...
                            m*obj.pixelScale*constants.radian2arcsec;
                        obj.frameHandle(1) = image(x,y,m_frame,...
                            'CDataMApping','Scaled',...
                            varargin{:});
                    end
                    obj.frameHandle(2) = title(sprintf('Frame #%d',obj.frameCount),...
                        'color',titleColor);
                    colormap(pink)
                    axis xy equal tight
                    colorbar('location','SouthOutside')
                    hu = findobj(gcf,'Type','uimenu','Label','OOMAO');
                    if isempty(hu)
                        hu = uimenu('Label','OOMAO');
                    end
                    hus  = uimenu(hu,'Label','Frame Listener Off','Callback',@oomaoMenu);
                    if obj.frameListener.Enabled
                        set(hus,'Label','Frame Listener On')
                    end
            end
            
            end
            function oomaoMenu(src,~)
                obj.frameListener.Enabled = ~obj.frameListener.Enabled;
                if obj.frameListener.Enabled
                    set(src,'Label','Frame Listener On')
                else
                    set(src,'Label','Frame Listener Off')
                end
            end
        end
        
        function varargout = grab(obj)
            %% GRAB Frame grabber
            %
            % grab(obj) grabs a frame
            %
            % out = grab(obj) grabs a frame and returns it
            
            switch class(obj.frameGrabber)
                case {'lensletArray','gpuLensletArray'}
                    readOut(obj,obj.frameGrabber.imagelets)
                case 'function_handle'
                    buffer = obj.frameGrabber();
                    [n,m] = size(buffer);
                    u = obj.roiSouthWestCorner(1):...
                        min(obj.roiSouthWestCorner(1)+obj.regionOfInterest(1)-1,n);
                    v = obj.roiSouthWestCorner(2):...
                        min(obj.roiSouthWestCorner(2)+obj.regionOfInterest(2)-1,m);
                    obj.frame = buffer(u,v);
                otherwise
            end
            if nargout>0
                varargout{1} = obj.frame;
            end
        end
        
        function relay(obj,src)
            
            % Here we check the last object the source went through before
            % the detector
            srcLastPath = src.opticalPath{end-1};
            switch class(srcLastPath) 
                case 'telescope'
                    f = utilities.cartAndPol(obj.resolution(1),...
                        'output','radius');
                    % pixel scale in radian
                    f = obj.pixelScale*f.*(obj.resolution(1)-1)./src.wavelength/2;
                    obj.frame = psf(srcLastPath,f);
                otherwise
                    if src.timeStamp>=obj.startDelay
                        obj.startDelay = -Inf;
                        obj.frameBuffer = obj.frameBuffer + src.intensity;
                        if src.timeStamp>=obj.exposureTime
                            src.timeStamp = 0;
                            disp(' @(detector:relay)> reading out and emptying buffer!')
                            readOut(obj,obj.frameBuffer)
                            obj.frameBuffer = 0*obj.frameBuffer;
                        end
                    end
            end
            
        end
        
    end
    
    methods (Access=protected)
        
        function readOut(obj,image)
            %% READOUT Detector readout
            %
            % readOut(obj,image) adds noise to the image: photon noise if
            % photonNoise property is true and readout noise if
            % readOutNoise property is greater than 0
            
            %             image = image;%This is now done in telescope.relay (.*obj.exposureTime;) % flux integration
            [n,m,~] = size(image);
            if any(n>obj.resolution(1))
%                 disp('Binning')
                image = utilities.binning(image,obj.resolution.*[n,m]/n);
            end
            obj.frameCount = obj.frameCount + 1;
            if obj.frameCount<obj.p_exposureTime*obj.clockRate
                obj.frameBuffer = obj.frameBuffer + image;
                obj.frame = [];
            else
                image = obj.frameBuffer + image;
                
                if license('checkout','statistics_toolbox')
                    if obj.photonNoise
                        image = poissrnd(image + obj.nPhotonBackground) - obj.nPhotonBackground;
                    end
                    image = obj.quantumEfficiency*image;
                    if obj.readOutNoise>0
                        image = image + randn(size(image)).*obj.readOutNoise;
                    end
                else
                    if obj.photonNoise
                        buffer    = image + obj.nPhotonBackground;
                        image = image + randn(size(image)).*(image + obj.nPhotonBackground);
                        index     = image<0;
                        image(index) = buffer(index);
                        image = obj.quantumEfficiency*image;
                    end
                    image = obj.quantumEfficiency*image;
                    if obj.readOutNoise>0
                        image = image + randn(size(image)).*obj.readOutNoise;
                    end
                end
                obj.frameCount = 0;
                obj.frame = image;
                obj.frameBuffer = 0;
            end
        end
        
    end
    
    methods (Static)
            
        function obj = loadobj(obj)
            %% LOADOBJ
            add(obj.log,obj,'Load!')
            setFrameListener(obj)
            obj.log = logBook.checkIn(obj);
        end
        
    end
    
    methods (Access=private)
        
        function setFrameListener(obj)
            %% SETFRAMELISTENER % Frame listener
            obj.frameListener = addlistener(obj,'frame','PostSet',...
                @(src,evnt) obj.imagesc );
            obj.frameListener.Enabled = false;
        end
        
    end
    
end
