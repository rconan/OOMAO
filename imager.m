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
        % Strehl ratio
        strehl;
        % entrapped energy
        ee;
        % entrapped energy slit width
        eeWidth;
    end
        
    properties (Access=private)
%         % integration count;
%         frameCount=0;
        % telescope
        tel;
    end
    
    methods
        
        %% Constructor
        function obj = imager(in,exposureTime,clockRate)
            if isa(in,'telescopeAbstract')
                resolution = in.resolution;
            elseif isnumeric(in)
                resolution = in;
            else
                error('oomao:imager','Inputer is either numeric or a telescope class')
            end
            obj = obj@detector(resolution);
            if isa(in,'telescopeAbstract')
                obj.tel = in;
                obj.exposureTime = in.samplingTime;
            end
            if nargin>1
                obj.exposureTime = exposureTime;
                obj.clockRate = clockRate;
            end
            obj.imgLens = lensletArray(1);
            obj.imgLens.nyquistSampling = 4;
            obj.imgLens.fieldStopSize = 10;
            % Frame listener
%             obj.frameListener = addlistener(obj,'frameBuffer','PostSet',...
%                 @(src,evnt) obj.imagesc );
%             obj.frameListener.Enabled = false;
        end
        
        function relay(obj,src)
            %% RELAY source propagation
            
            relay(obj.imgLens,src)
%             if all([src.timeStamp]>=obj.startDelay)
%                 if obj.frameCount==0
%                     disp(' @(detector:relay)> starting imager integration!')
%                 end
%                 obj.startDelay = -Inf;
%                 obj.frameCount = obj.frameCount + 1;
%                 obj.frameBuffer = obj.frameBuffer + cat(2,src.intensity);
%                 if src.timeStamp>=obj.exposureTime
%                     src.timeStamp = 0;
%                     flush(obj,src);
%                 end
%             end
              readOut(obj,obj.imgLens.imagelets)
              if obj.frameCount==0 && obj.startDelay==0
                  flush(obj)
              end
        end
        
        function flush(obj)
            fprintf(' @(detector:relay)> reading out and emptying buffer (%d frames)!\n',obj.frameCount)
            if ~isempty(obj.referenceFrame) && ~isempty(obj.frame)
                obj.imgLens.fieldStopSize = obj.imgLens.fieldStopSize*2;
                src_ = source.*obj.referenceFrame;
                wavePrgted = propagateThrough(obj.imgLens,src_);
                otf =  abs(wavePrgted);
                nFrame = obj.exposureTime*obj.clockRate;
                
                n = length(obj.referenceFrame);
                m_frame = obj.frame/nFrame;
                nf = size(m_frame)/n;                
                m_frame = mat2cell( m_frame , n*ones(1,nf(1)), n*ones(1,nf(2)));
                
                obj.strehl = zeros(1,length(m_frame));
                obj.ee     = zeros(1,length(m_frame));
                a = [];
                for kFrame = 1:length(m_frame)
                    src_ = src_.*m_frame{kFrame};
                    wavePrgted = propagateThrough(obj.imgLens,src_);
                    otfAO =  abs(wavePrgted);
                    a = [a otfAO];
                    % strehl ratio
                    obj.strehl(kFrame) = sum(otfAO(:))/sum(otf(:));
                    % entrapped energy
%                     if ~isempty(obj.eeWidth)
%                         a      = (obj.eeWidth/(src_.wavelength/obj.tel.D*constants.radian2arcsec))/obj.tel.D;
%                         nOtf   = length(otfAO);
%                         u      = linspace(-1,1,nOtf).*obj.tel.D;
%                         [x,y]  = meshgrid(u);
%                         eeFilter ...
%                             = a^2*(sin(pi.*x.*a)./(pi.*x.*a)).*...
%                             (sin(pi.*y.*a)./(pi.*y.*a));
%                         otfAO = otfAO/max(otfAO(:));
%                         obj.ee(kFrame) = real(trapz(u,trapz(u,otfAO.*eeFilter)));
%                     end

                end
                obj.imgLens.fieldStopSize = obj.imgLens.fieldStopSize/2;
            end
            obj.frameCount = 0;
        end
        
%         function imagesc(obj,varargin)
%             %% IMAGESC Display the detector frame
%             %
%             % imagesc(obj) displays the frame of the detector object
%             %
%             % imagesc(obj,'PropertyName',PropertyValue) displays the frame of
%             % the detector object and set the properties of the graphics object
%             % imagesc
%             %
%             % h = imagesc(obj,...) returns the graphics handle
%             %
%             % See also: imagesc
%             
%             disp('Got it')
%             if ishandle(obj.frameHandle)
%                 set(obj.frameHandle,'Cdata',obj.frameBuffer,varargin{:});
%                 %                 xAxisLim = [0,size(obj.frame,2)]+0.5;
%                 %                 yAxisLim = [0,size(obj.frame,1)]+0.5;
%                 %                 set( get(obj.frameHandle,'parent') , ...
%                 %                     'xlim',xAxisLim,'ylim',yAxisLim);
%             else
%                 obj.frameHandle = image(obj.frameBuffer,...
%                     'CDataMApping','Scaled',...
%                     varargin{:});
%                 colormap(pink)
%                 axis xy equal tight
%                 colorbar('location','SouthOutside')
%             end
%         end

   
    end

end