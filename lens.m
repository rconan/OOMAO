classdef lens < handle
% LENS Create a lens array object
%
% obj = lensletArray(nLenslet) creates a lenslet array object from the
% number of lenslet on one side of the array

    properties
        % the # of lenslet on one side of the lenslet array
        nLenslet;
        % lenslet size
        pitch;
        % the minimum amount of light per lenslet
        minLightRatio;
        % lenslets conjugation altitude
        conjugationAltitude = Inf;
        % the input wave
        wave;
        % imagelets listener
        imageletsListener
        % stacked imagelets sum
        sumStack = false;
        % the number of lenslet arrays
        nArray;
        % optical throughput
        throughput=1;
        % lenslet array tag
        tag='LENSLET ARRAY';
    end
    
    properties (SetObservable=true)
        % the image at each lenslet focus
        imagelets;
    end
    
    properties (SetAccess=private)
        fftPad;        
    end

    properties (Dependent)
        % the # of pixel per lenslet
        nLensletWavePx;
        % the Nyquist sampling factor
        nyquistSampling;
        % the lenslet field of view given in diffraction fwhm units
        fieldStopSize;
    end
    
    properties (Dependent, SetAccess=private)
        % the # of pixel per imagelet
        nLensletImagePx;
        % the total # of pixel per imagelet
        nLensletsImagePx;
    end

    properties (Access=private)
        p_nLensletWavePx;
        p_nyquistSampling;
        p_fieldStopSize;
        fftPhasor;
        imageHandle;
        zoomTransmittance;
        log;
    end

    methods

        %% Constructor
        function obj = lens(varargin)
%             error(nargchk(1, 3, nargin))
            obj.nLenslet        = 1;
            obj.minLightRatio   = 0;
            obj.nyquistSampling = 1;
            obj.imageletsListener = addlistener(obj,'imagelets','PostSet',...
                @(src,evnt) obj.imagesc );
            obj.imageletsListener.Enabled = false;
            obj.log = logBook.checkIn(obj);
        end

        %% Destructor
        function delete(obj)
            if ishandle(obj.imageHandle)
                delete(get(obj.imageHandle,'parent'));
            end
            checkOut(obj.log,obj)
        end
        
        function display(obj)
            %% DISPLAY Display object information
            %
            % disp(obj) prints information about the lenslet array object
          
            fprintf('___ %s ___\n',obj.tag)
            if obj.nArray>1
                fprintf(' %d %dx%d lenslet array: \n',...
                    obj.nArray,obj.nLenslet*ones(1,2))
            else
                fprintf(' %dx%d lenslet array: \n',...
                    obj.nLenslet*ones(1,2))
            end
            fprintf('  . %3.1f pixels across the diffraction limited spot fwhm\n',...
                obj.nyquistSampling*2);
            fprintf('  . %d pixels across the square lenslet field stop size\n',...
                obj.fieldStopSize*obj.nyquistSampling*2);            
            fprintf('  . optical throughput coefficient: %3.1f\n',...
                obj.throughput);            
            fprintf('----------------------------------------------------\n')
            
        end

        %% Lenslet image sampling
        function nLensletImagePx = get.nLensletImagePx(obj)
            nLensletImagePx = ceil(obj.fieldStopSize.*obj.nyquistSampling*2);
        end

        %% Whole lenslet array image sampling
        function nLensletsImagePx = get.nLensletsImagePx(obj)
            nLensletsImagePx = obj.nLenslet.*obj.nLensletImagePx;
        end

        %% Nyquist sampling
        function out = get.nyquistSampling(obj)
            out = obj.p_nyquistSampling;
        end
        function set.nyquistSampling(obj,val)
            obj.p_nyquistSampling = val;
            obj.fftPad = ceil(obj.p_nyquistSampling*2);
            if ~isempty(obj.nLensletWavePx) % refreshing phasor
                setPhasor(obj);
            end
        end
        
        %% Field stop
        function out = get.fieldStopSize(obj)
            out = obj.p_fieldStopSize;
        end
        function set.fieldStopSize(obj,val)
            obj.p_fieldStopSize = val;
            if ~isempty(obj.nLensletWavePx) % refreshing phasor
                setPhasor(obj);
            end
        end

        %% Number of wave pixel per lenslet
        function out = get.nLensletWavePx(obj)
            out = obj.p_nLensletWavePx;
        end
        function set.nLensletWavePx(obj,val)
            obj.p_nLensletWavePx = val;
                fprintf(' @(lensletArray)> Setting the lenslet field stop size!\n')
                obj.fieldStopSize  = obj.p_nLensletWavePx/obj.p_nyquistSampling/2;
%             % set the wave reshaping index (2D to 3D)
%             logBook.add(obj,'Set the wave reshaping index (2D to 3D)')
        end  

        function relay(obj,src)
            %% PROPAGATETHROUGH Fraunhoffer wave propagation to the lenslets focal plane
            %
            % propagateThrough(obj) progates the object wave throught the
            % lenslet array and computes the imagelets
            
            % for a given wave of size [nxn], 
            % nyquistSampling=1 means 2 pixels per fhwm obtained by setting
            % fftPad=2, fieldStopSize should be set to default n/nLenslet
            % fwhm
            % nyquistSampling=0.5 means 1 pixels per fhwm obtained by
            % setting fftPad=1, fieldStopSize should be set to default
            % n/nLenslet/nyquistSampling/2
            
            if ndims(src)==3
                obj.sumStack = true;
                obj.nArray = 1;
            else
                obj.nArray = numel(src);
            end
            val = src.catWave;
            if isscalar(val) % if source amplitude and phase not set, set a default one
                n = obj.nLensletWavePx*obj.nLenslet;
                set(src,...
                    'mask',true(n),...
                    'amplitude',ones(n),...
                    'phase',zeros(n));
                if ~isempty(src(1).opticalPath)
                    cellfun(@(x)relay(x,src),src(1).opticalPath,'uniformOutput',false)
                end
                val = src.catWave;
            end
            [nLensletsWavePx,nLensletsWavePxNGuideStar,nWave] = size(val);
            % Resize the 3D input into a 2D input
            nLensletsWavePxNGuideStar = nLensletsWavePxNGuideStar*nWave;
            val = reshape(val,nLensletsWavePx,nLensletsWavePxNGuideStar);
            nLensletWavePx = nLensletsWavePx./obj.nLenslet;
            if isempty(obj.nLensletWavePx) || all(obj.nLensletWavePx~=nLensletWavePx)
                obj.nLensletWavePx = nLensletWavePx;
            end
            nLensletArray = nLensletsWavePxNGuideStar/nLensletsWavePx;
            % Invocation of the zoom optics for conjugation to finite
            % distance
%             if isfinite(obj.conjugationAltitude)
% %                 val                 = obj.zoomTransmittance.*val;
%                 val = repmat(obj.zoomTransmittance,1,nLensletArray).*val;
%             end
%             nOutWavePx    = obj.nLensletImagePx*obj.fftPad;    % Pixel length of the output wave
%             evenOdd       = rem(obj.nLensletImagePx,2);
            nOutWavePx    = obj.nLensletWavePx*obj.fftPad;    % Pixel length of the output wave
            evenOdd       = rem(obj.nLensletWavePx,2);
            if ~rem(nOutWavePx,2) && evenOdd
                nOutWavePx = nOutWavePx + evenOdd;
            end
            nOutWavePx = max(nOutWavePx,obj.nLensletWavePx);
            nLensletSquareWavePx    = obj.nLensletWavePx*obj.nLenslet^2*nLensletArray;
            wavePrgted = zeros(nOutWavePx,nLensletSquareWavePx);
            val        = val./nOutWavePx;
%             nLensletWavePx   = obj.nLensletWavePx;
            nLensletImagePx  = obj.nLensletImagePx;
            nLensletsImagePx = obj.nLensletsImagePx;
            %%% ODD # OF PIXELS PER LENSLET
            if isempty(obj.fftPhasor)
%                 disp('Phasor empty!')
                % Shape the wave per columns of lenslet pixels
                val       = reshape(val,obj.nLensletWavePx,nLensletSquareWavePx);
                u         = any(val); % Index of non-zeros columns
                wavePrgted(:,u) = fftshift(fft(val(:,u),nOutWavePx),1);
                % Select the field of view
                v = [];
                if nOutWavePx>nLensletImagePx
%                     centerIndex = (nOutWavePx+1)/2;
%                     halfLength  = (nLensletImagePx-1)/2;
                    centerIndex = ceil((nOutWavePx+1)/2);
                    halfLength  = floor(nLensletImagePx/2);
                    v           = true(nOutWavePx,1);
                    v((-halfLength:halfLength)+centerIndex) ...
                                = false;
                elseif nOutWavePx<nLensletImagePx
                    error('lensletArray:propagateThrough:size','The computed image is smaller than the expected image!')
                end
                wavePrgted(v,:) = [];
                % Back to transpose 2D
                val       = reshape( wavePrgted ,...
                    nLensletsImagePx,obj.nLensletWavePx*obj.nLenslet*nLensletArray).';
                % Shape the wave per rows of lenslet pixels
                val       = reshape(val,obj.nLensletWavePx,nLensletsImagePx*obj.nLenslet*nLensletArray);
                u         = any(val); % Index of non-zeros columns
                wavePrgted = zeros(nOutWavePx,nLensletsImagePx*obj.nLenslet*nLensletArray);
                wavePrgted(:,u)  = fftshift(fft(val(:,u),nOutWavePx),1);
                wavePrgted(v,:) = [];
            else
            %%% EVEN # OF PIXELS PER LENSLET
%                 disp('Phasor exist!')
                val       = val.*repmat(obj.fftPhasor,1,nLensletArray);
                % Shape the wave per columns of lenslet pixels
                val       = reshape(val,obj.nLensletWavePx,nLensletSquareWavePx);
                u         = any(val); % Index of non-zeros columns
                wavePrgted(:,u) = fft(val(:,u),nOutWavePx);
                v = [];
                if nOutWavePx>nLensletImagePx
%                     disp('Cutting through!')
%                     centerIndex = nOutWavePx/2+1;
%                     halfLength  = nLensletImagePx/2;
                    centerIndex = ceil((nOutWavePx+1)/2);
                    halfLength  = floor(nLensletImagePx/2);
                    v           = true(nOutWavePx,1);
                    v((0:nLensletImagePx-1)+centerIndex-halfLength) ...
                                = false;
                elseif nOutWavePx<nLensletImagePx;
                    error('lensletArray:propagateThrough:size','The computed image is smaller than the expected image!')
                end
                wavePrgted(v,:) = [];
                % Back to transpose 2D
                val       = reshape( wavePrgted ,...
                    nLensletsImagePx,obj.nLensletWavePx*obj.nLenslet*nLensletArray).';
                % Shape the wave per rows of lenslet pixels
                val       = reshape(val,obj.nLensletWavePx,nLensletsImagePx*obj.nLenslet*nLensletArray);
                u         = any(val); % Index of non-zeros columns
                wavePrgted = zeros(nOutWavePx,nLensletsImagePx*obj.nLenslet*nLensletArray);
                wavePrgted(:,u)  = fft(val(:,u),nOutWavePx);
                wavePrgted(v,:) = [];
            end
            % Back to transpose 2D
            wavePrgted  = reshape(wavePrgted,nLensletsImagePx*nLensletArray,nLensletsImagePx).';
%             wavePrgted = wavePrgted.*conj(wavePrgted);
            % and back to input wave array shape
            [n,m] = size(wavePrgted);
            wavePrgted = reshape(wavePrgted,n,n,[]);
%             if obj.sumStack
%                 wavePrgted = sum(wavePrgted,3);
%             end
            src = resetAmplitudeAndPhase(src);
            for kSrc = 1:length(src)
            src(kSrc).amplitude = abs(wavePrgted(:,:,kSrc));
            src(kSrc).phase     = angle(wavePrgted(:,:,kSrc));
            src(kSrc).mask = true(n);
            end
            obj.imagelets = wavePrgted.*conj(wavePrgted)*obj.throughput;
        end
        
%         function relay(obj,src)
%             propagateThrough(obj,src)
%         end

        function varargout = imagesc(obj,varargin)
        %% IMAGESC Display the lenslet imagelets
        %
        % imagesc(obj) displays the imagelets of the lenslet array object 
        %
        % imagesc(obj,'PropertyName',PropertyValue) displays the imagelets
        % of the lenslet array object and set the properties of the
        % graphics object imagesc
        %
        % h = imagesc(obj,...) returns the graphics handle
        %
        % See also: imagesc
        
 %             if ishandle(obj.imageHandle)
%                 axisLim = [0,obj.nLensletsImagePx]+0.5;
%                 set( get(obj.imageHandle,'parent') , ...
%                     'xlim',axisLim,'ylim',axisLim);
%             end
       if ishandle(obj.imageHandle)
            set(obj.imageHandle,'Cdata',obj.imagelets,varargin{:});
        else
            %                 figure
            obj.imageHandle = image(obj.imagelets,...
                'CDataMApping','Scaled',varargin{:});
            colormap(pink)
            axis xy square
            axisXLim = [0,size(obj.imagelets,2)]+0.5;
            axisYLim = [0,size(obj.imagelets,1)]+0.5;
            set( get(obj.imageHandle,'parent') , ...
                'xlim',axisXLim,'ylim',axisYLim);
            axis equal tight
            colorbar
        end
        if nargout>0
            varargout{1} = obj.imageHandle;
        end
        end


    end
    
    methods (Access=private)     
        
        function setPhasor(obj)
        % Phasor settings for Fraunhoffer propagation in order to have the
        % spots centered in between pixels for even pixel sampling
%             nOutWavePx    = obj.nLensletImagePx*obj.fftPad;    % Pixel length of the output wave
%             evenOdd       = rem(obj.nLensletImagePx,2);
%             if ~rem(nOutWavePx,2) && evenOdd
%                 nOutWavePx = nOutWavePx + evenOdd;
%             end
%             fprintf(' 0ld nOutWavePx = %d',nOutWavePx);
            nOutWavePx    = obj.nLensletWavePx*obj.fftPad;    % Pixel length of the output wave
            evenOdd       = rem(obj.nLensletWavePx,2);
            if ~rem(nOutWavePx,2) && evenOdd
                nOutWavePx = nOutWavePx + evenOdd;
            end
            nOutWavePx = max(nOutWavePx,obj.nLensletWavePx);
%             fprintf(' - new nOutWavePx = %d\n',nOutWavePx);
            if ~rem(obj.nLensletImagePx,2)
                % shift the intensity of half a pixel for even sampling
                fprintf(' @(lensletArray)> Set phasor (shift the intensity of half a pixel\n for even intensity sampling)\n')
                [u,v]         = ndgrid((0:(obj.nLensletWavePx-1)).*(1-nOutWavePx)./nOutWavePx);
                obj.fftPhasor = repmat( exp(-1i.*pi.*(u+v)) , obj.nLenslet, obj.nLenslet );
            else
                fprintf(' @(lensletArray)> Reset phasor\n')
                obj.fftPhasor = [];
            end            
        end

    end

end