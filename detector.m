classdef detector < handle
    % DETECTOR Create a detector object
    %
    % obj = detector(resolution) creates a detector object from the detector
    % resolution
    
    properties
        % Add or not photon noise to the frame
        photonNoise = false;
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
        % Exposure time [second]
        exposureTime = 1;
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
        frameBuffer = 0;;
    end
    
    properties (SetObservable=true)
        % detector frame
        frame;
    end
    
    properties (Access=protected)
        frameHandle;
        log;
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
            % Frame listener
            obj.frameListener = addlistener(obj,'frame','PostSet',...
                @(src,evnt) obj.imagesc );
            obj.frameListener.Enabled = false;
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
            display(obj)
        end
        
        %% Destructor
        function delete(obj)
            if ishandle(obj.frameHandle)
                delete(get(obj.frameHandle,'Parent'));
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
            fprintf('  . %.1f photo-events rms read-out moise \n',...
                obj.readOutNoise)
            fprintf('  . %3.1fms exposure time and %3.1fHz frame rate \n',...
                obj.exposureTime*1e3,obj.frameRate)
            fprintf('----------------------------------------------------\n')
            
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
            
            if ishandle(obj.frameHandle)
                set(obj.frameHandle,'Cdata',obj.frame,varargin{:});
                %                 xAxisLim = [0,size(obj.frame,2)]+0.5;
                %                 yAxisLim = [0,size(obj.frame,1)]+0.5;
                %                 set( get(obj.frameHandle,'parent') , ...
                %                     'xlim',xAxisLim,'ylim',yAxisLim);
            else
                if isempty(obj.pixelScale)
                    obj.frameHandle = image(obj.frame,...
                        'CDataMApping','Scaled',...
                        varargin{:});
                else
                    [n,m] = size(obj.frame);
                    x = 0.5*linspace(-1,1,n)*...
                        n*obj.pixelScale*constants.radian2arcsec;
                    y = 0.5*linspace(-1,1,m)*...
                        m*obj.pixelScale*constants.radian2arcsec;
                    obj.frameHandle = image(x,y,obj.frame,...
                        'CDataMApping','Scaled',...
                        varargin{:});
                end
                colormap(pink)
                axis xy equal tight
                colorbar('location','SouthOutside')
            end
        end
        
        function varargout = grab(obj)
            %% GRAB Frame grabber
            %
            % grab(obj) grabs a frame
            %
            % out = grab(obj) grabs a frame and returns it
            
            switch class(obj.frameGrabber)
                case 'lensletArray'
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
            [n,m,k] = size(image);
            if any(n>obj.resolution(1))                
%                 disp('Binning')
                image = utilities.binning(image,obj.resolution);
            end
            
            if license('checkout','statistics_toolbox')
                if obj.photonNoise
                    image = poissrnd(image);
                end
                image = obj.quantumEfficiency*image;
                if obj.readOutNoise>0
                    image = image + randn(size(image)).*obj.readOutNoise;
                end
            else
                if obj.photonNoise
                    buffer    = image;
                    image = image + randn(size(image)).*image;
                    index     = image<0;
                    image(index) = buffer(index);
                    image = obj.quantumEfficiency*image;
                end
                image = obj.quantumEfficiency*image;
                if obj.readOutNoise>0
                    image = image + randn(size(image)).*obj.readOutNoise;
                end
            end
            obj.frame = image;
        end
        
    end
    
end
