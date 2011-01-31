classdef deformableMirror < handle
    % DEFORMABLEMIRROR Create a deformableMirror object
    %
    % dm = deformableMirror(nActuator,'modes',ifObject) creates a
    % deformableMirror object from the number of actuator across the mirror
    % diameter and from the influence function object
    %
    % dm = deformableMirror(nActuator,'modes',IF) creates a
    % deformableMirror object from the number of actuator across the mirror
    % diameter and from the influence function matrix
    %
    % See also influenceFunction
   
    properties
        % conjugation altitude
        zLocation = 0;
        % # of actuator
        nActuator;
        % influence functions
        modes;
        % electronic driver
        driver;
        % coefficients default
        coefsDefault;
        % surface listener
        surfaceListener;
        % lexicographic ordering of DM surface map (default: false)
        lex = false;
        % units of dm coefficients [default: micron]
        coefsUnit = 1;%e-6;
        % deformableMirror tag
        tag = 'DEFORMABLE MIRROR';
    end
    
    properties (SetObservable=true,Dependent,SetAccess=private)
        % the shape of the DM
        surface;
    end
    
    properties (Dependent,SetAccess=private)
         % # of actuators in the pupil
        nValidActuator;
   end
    
    properties (Dependent)
        % valid actuator
        validActuator;
   end
    
    properties (Dependent, SetObservable=true)
        % coefficients
        coefs;
    end
    
    properties (Access=private)
        p_coefs; 
        p_surface;
        p_validActuator;
        imageHandle;
        log;
    end
        
    methods
        
        %% Constructor
        function obj = deformableMirror(nActuator,varargin)
            p = inputParser;
            p.addRequired('nActuator', @isnumeric);
            p.addParamValue('modes', [], @(x) isnumeric(x) || (isa(x,'influenceFunction') || isa(x,'zernike')) );
            p.addParamValue('resolution', [], @isnumeric);
            p.addParamValue('validActuator', ones(nActuator), @islogical);
            p.addParamValue('zLocation', 0, @isnumeric);
            p.parse(nActuator, varargin{:});
            obj.nActuator         = p.Results.nActuator;
            obj.p_validActuator     = p.Results.validActuator;
            obj.modes             = p.Results.modes;
            obj.zLocation             = p.Results.zLocation;
            obj.surfaceListener = addlistener(obj,'surface','PostSet',...
                @(src,evnt) obj.imagesc );
            obj.surfaceListener.Enabled = false;
            if isa(obj.modes,'influenceFunction') && ~isempty(p.Results.resolution)
                setInfluenceFunction(obj.modes,obj.nActuator,p.Results.resolution,obj.validActuator,1,[0,0]);
            elseif isa(obj.modes,'zernike')
                obj.p_validActuator = true(1,obj.modes.nMode);
            end
            obj.coefsDefault      = zeros(obj.nValidActuator,1);
            obj.coefs             = zeros(obj.nValidActuator,1);
            obj.log = logBook.checkIn(obj);
            display(obj)
        end
        
        %% Destructor
        function delete(obj)
%             if isa(obj.modes,'influenceFunction')
%                 delete(obj.modes)
%             end
            checkOut(obj.log,obj)
        end
        
        function display(obj)
            %% DISPLAY Display object information
            %
            % display(obj) prints information about the deformable mirror
            % object
          
            fprintf('___ %s ___\n',obj.tag)
            fprintf(' %dX%d actuators deformable mirror: \n  . %d controlled actuators\n',...
                obj.nActuator,obj.nActuator,obj.nValidActuator)
            fprintf('----------------------------------------------------\n')
            if isa(obj.modes,'influenceFunction')
                display(obj.modes)
            end

        end
        
        %% Get nValidActuator
        function out = get.nValidActuator(obj)
            out = sum(obj.validActuator(:));
        end
        
        %% Set and Get coefs
        function out = get.coefs(obj)
            out = obj.p_coefs;
        end
        function set.coefs(obj,val)
            if isscalar(val)
                val = ones(obj.nValidActuator,1)*val;
            end
            obj.p_coefs = obj.coefsUnit*bsxfun(@plus,val,obj.coefsDefault);
            if isa(obj.driver,'function_handle')
                obj.driver(obj.p_coefs);
                return
            end
            if isvalid(obj.surfaceListener) && obj.surfaceListener.Enabled
                imagesc(obj);
            end
        end
        
        %% Set and Get validActuator
        function out = get.validActuator(obj)
            out = obj.p_validActuator;
        end
        function set.validActuator(obj,val)
            if obj.nValidActuator>=sum(val(:))
                val(~obj.p_validActuator)  = [];
                obj.coefsDefault(~val)     = [];
                obj.p_coefs(~val)          = [];
                obj.modes.modes(:,~val)    = [];
                obj.p_validActuator        = val;
            else
                error('oomao:deformableMirror:set.validActuator','Sorry only downsizing is possible!')
            end
        end        
        
        %% Get the dm shape
        function dmShape = get.surface(obj)
            obj.p_surface = obj.modes*obj.coefs;
            if obj.lex
                dmShape = obj.p_surface;
            else
                dmShape = utilities.toggleFrame(obj.p_surface,3);
            end
        end
        
        function relay(obj,src)
            %% RELAY deformable mirror to source relay
            %
            % relay(obj,srcs) writes the deformableMirror amplitude and
            % phase into the properties of the source object(s)
            
            nSrc       = numel(src);
            wavenumber = 2*pi/src(1).wavelength;
            dmPhase    = -2*obj.surface*wavenumber;
            nPhase     = size(obj.p_coefs,2);
            if nPhase>nSrc
                for kSrc = 1:nSrc
                    src(kSrc).phase = dmPhase;
                    src(kSrc).amplitude = 1;
                end
            else
                for kSrc = 1:nSrc
                    src(kSrc).phase = dmPhase(:,:,min(kSrc,nPhase));
                    src(kSrc).amplitude = 1;
                end
            end
        end
        
        function out = fittingError(obj,telAtm,src,unit)
            %% FITTINGERROR deformable mirror fitting error
            %
            % out = fittingError(telAtm) computes the deformable mirror
            % fitting error variance in radian^2 for the given
            % telescope+atmosphere system 
            % out = fittingError(telAtm,src) computes the deformable mirror
            % fitting error rms in meter for the given telescope+atmosphere
            % system
            % out = fittingError(telAtm,src,unit) computes the deformable
            % mirror fitting error rms in meterX10^-unit for the given
            % telescope+atmosphere system (nanometer: unit=-9)
            
            add(obj.log,obj,' Computing the fitting error ...')
            atm = telAtm.opticalAberration;
            atmWavelength  = atm.wavelength;
            atm.wavelength = src.wavelength; 
            if isa(obj.modes,'zernike')
                out = zernikeStats.residualVariance(obj.modes.nMode,atm,telAtm);
            else
                d = telAtm.D/(obj.nActuator-1);
                fc = 1/d/2;
                a = phaseStats.variance(telAtm.opticalAberration);
                b = dblquad( @(fx,fy) phaseStats.spectrum( hypot(fx,fy) , atm ) , ...
                    -fc,fc,-fc,fc);
                out = a - b;
            end
            atm.wavelength = atmWavelength;
            if nargin>2
                if nargin<4
                    unit = 1;
                end
                out = 10^-unit*sqrt(out)/src.waveNumber;
            end            
        end
        
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
            
            if ishandle(obj.imageHandle)
                set(obj.imageHandle,'Cdata',obj.surface,varargin{:});
            else
                %                 figure
                obj.imageHandle = image(obj.surface,...
                    'CDataMApping','Scaled',varargin{:});
                set( get(obj.imageHandle,'parent') ,'visible','off')
                colormap(pink)
                axis square
                colorbar
            end
            if nargout>0
                varargout{1} = obj.imageHandle;
            end
        end
        
        
    end
    
end