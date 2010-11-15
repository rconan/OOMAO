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
        % entrapped energy
        ee;
        % entrapped energy slit width
        eeWidth;
    end
        
    properties (Access=private)
        % integration count;
        frameCount=0;
        % telescope
        tel;
    end
    
    methods
        
        %% Constructor
        function obj = imager(in)
            if isa(in,'telescope')
                resolution = in.resolution;
            elseif isnumeric(in)
                resolution = in;
            else
                error('oomao:imager','Inputer is either numeric or a telescope class')
            end
            obj = obj@detector(resolution);
            if isa(in,'telescope')
                obj.tel = in;
                obj.exposureTime = in.samplingTime;
            end
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
                        obj.imgLens.fieldStopSize = obj.imgLens.fieldStopSize*2;
                        src_ = source.*obj.referenceFrame*obj.imgLens;
                        otf =  src_.amplitude;
                        src_ = src_.*(obj.frame/obj.frameCount)*obj.imgLens;
                        obj.imgLens.fieldStopSize = obj.imgLens.fieldStopSize/2;
                        otfAO =  src_.amplitude;
%                         figure, imagesc(real(otfAO)/max(otfAO(:)))
                        % strehl ratio
                        obj.strehl = sum(otfAO(:))/sum(otf(:));
                        % entrapped energy
                        a      = (obj.eeWidth/(src.wavelength/obj.tel.D*constants.radian2arcsec))/obj.tel.D;
                        nOtf   = length(otfAO);
                        u      = linspace(-1,1,nOtf).*obj.tel.D;
                        [x,y]  = meshgrid(u);
                        eeFilter ...
                               = a^2*(sin(pi.*x.*a)./(pi.*x.*a)).*...
                            (sin(pi.*y.*a)./(pi.*y.*a));
                        otfAO = otfAO/max(otfAO(:));
                        obj.ee = real(trapz(u,trapz(u,otfAO.*eeFilter)));
                    end
                end
                disp(obj.frameCount)
                obj.frameCount = 0;
            end
        end
   
    end

end