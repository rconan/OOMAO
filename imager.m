classdef imager < detector
    %% Imaging camera
    %
    % imgr = imager(resolution) creates an imaging camera object from the
    % detector resolution
    %
    % Example:
    % tel = telescope(1,'resolution',21,'samplingTime',1);
    % imgr = imager(21);
    % src = source.*tel*imgr;
    % figure
    % imagesc(imgr.frame)

    properties
        % reference frame
        referenceFrame;
        % imaging lens
        imgLens;
        % imaging camera
        imgCcd;
        % Strehl ratio
        strehl;
    end
    
    properties (Access=private)
        % integration count;
        frameCount=0;
    end
    
    methods
        
        %% Constructor
        function obj = imager(resolution)
            obj = obj@detector(resolution);
            obj.imgLens = lens;
        end
        
        function relay(obj,src)
            %% RELAY source propagation
            
            relay(obj.imgLens,src)
            if src.timeStamp>=obj.startDelay
                obj.startDelay = -Inf;
                obj.frameCount = obj.frameCount + 1;
                obj.frameBuffer = obj.frameBuffer + src.intensity;
                if src.timeStamp>=obj.exposureTime
                    src.timeStamp = 0;
                    disp(' @(detector:relay)> reading out and emptying buffer!')
                    readOut(obj,obj.frameBuffer)
                    obj.frameBuffer = 0*obj.frameBuffer;
                    if ~isempty(obj.referenceFrame)
                        src_ = source.*obj.referenceFrame*obj.imgLens;
                        otf =  src_.amplitude;
                        src_ = src_.*(obj.frame/obj.frameCount)*obj.imgLens;
                        otfAO =  src_.amplitude;
                        obj.strehl = sum(otfAO(:))/sum(otf(:));
                    end
                end
                obj.frameCount = 0;
            end
        end
   
    end

end