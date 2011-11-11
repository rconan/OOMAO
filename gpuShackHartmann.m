classdef gpuShackHartmann < shackHartmann
    
    methods
        
        %% Constructor
        function obj = gpuShackHartmann(nLenslet,detectorResolution,minLightRatio)
            obj = obj@shackHartmann;
            error(nargchk(1, 4, nargin))
            obj.lenslets = gpuLensletArray(nLenslet);
            obj.camera   = detector(detectorResolution);
            if detectorResolution==2
                obj.quadCell = true;
                obj.centroiding = false;
            end
            obj.lenslets.nLensletWavePx = ...
                detectorResolution/nLenslet;
            if nargin>2
                obj.lenslets.minLightRatio = minLightRatio;
            else
                obj.lenslets.minLightRatio = 0;
            end
            obj.validLenslet = true(nLenslet);
            obj.camera.frameGrabber ...
                = obj.lenslets;
            obj.referenceSlopes = zeros(obj.nValidLenslet*2,1);
            obj.p_referenceSlopes = ...
                repmat(obj.p_referenceSlopes,obj.lenslets.nArray,1);
            
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
            display(obj)
        end
        
    end
    
end