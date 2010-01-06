classdef detector < handle
    % DETECTOR Create a detector object
    %
    % obj = detector(resolution) creates a detector object from the detector
    % resolution
    
    properties
        % Add or not photon noise to the frame
        photonNoiseLess = true;
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
        % detector timer
        paceMaker;
        % frame update listener
        frameListener;
        % frame grabber callback function
        frameGrabber;
    end
    
    properties (SetObservable=true)
        % detector frame
        frame;
    end
    
    properties (GetAccess=private,SetAccess=private)
        frameHandle;
    end
    
    methods
        
        % Constructor
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
            obj.paceMaker.TimerFcn = {@timerCallBack, obj};
            obj.paceMaker.ExecutionMode = 'FixedSpacing';
            %             obj.paceMaker.BusyMode = 'error';
            obj.paceMaker.Period = 1;
            obj.paceMaker.ErrorFcn = 'disp('' @detector: frame rate too high!'')';
            function timerCallBack( timerObj, event, a)
                %                 fprintf(' @detector: %3.2fs\n',timerObj.instantPeriod)
                a.grab;
            end
            %             obj.frameRate = 1;
        end
        
        % Destructor
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
        end
        
        function imagesc(obj,varargin)
            % IMAGESC Display the detector frame
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
                obj.frameHandle = image(obj.frame,...
                    'CDataMApping','Scaled',...
                    varargin{:});
                colormap(pink)
                axis equal tight
                colorbar
            end
        end
        
        function varargout = grab(obj)
            % GRAB Frame grabber
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
        
    end
    
    methods (Access=protected)
        
        function readOut(obj,image)
            % READOUT Detector readout
            %
            % readOut(obj,image) adds noise to the image: photon noise if
            % photonNoiseLess property is false and readout noise if
            % readOutNoise property is greater than 0
            
            image = image.*obj.exposureTime; % flux integration
            if license('checkout','statistics_toolbox')
                if ~obj.photonNoiseLess
                    image = obj.quantumEfficiency*poissrnd(image);
                end
                if obj.readOutNoise>0
                    image = normrnd(image,obj.readOutNoise);
                end
            else
                if ~obj.photonNoiseLess
                    buffer    = image;
                    image = image + randn(size(image)).*image;
                    index     = image<0;
                    image(index) = buffer(index);
                    image = obj.quantumEfficiency*image;
                end
                if obj.readOutNoise>0
                    image = image + randn(size(image)).*obj.readOutNoise;
                end
            end
            obj.frame = image;
        end
        
    end
    
end
