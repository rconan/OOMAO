classdef deformableMirror < handle
   
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
        % valid actuator
        validActuator;
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
    
    properties (Dependent, SetObservable=true)
        % coefficients
        coefs;
    end
    
    properties (Access=private)
        p_coefs; 
        p_surface;
        imageHandle;
        log;
    end
        
    methods
        
        %% Constructor
        function obj = deformableMirror(nActuator,varargin)
            p = inputParser;
            p.addRequired('nActuator', @isnumeric);
            p.addParamValue('modes', [], @(x) isnumeric(x) || isa(x,'influenceFunction') );
            p.addParamValue('resolution', [], @isnumeric);
            p.addParamValue('validActuator', ones(nActuator), @islogical);
            p.addParamValue('zLocation', 0, @isnumeric);
            p.parse(nActuator, varargin{:});
            obj.nActuator         = p.Results.nActuator;
            obj.validActuator     = p.Results.validActuator;
            obj.modes             = p.Results.modes;
            obj.zLocation             = p.Results.zLocation;
            obj.surfaceListener = addlistener(obj,'surface','PostSet',...
                @(src,evnt) obj.imagesc );
            obj.surfaceListener.Enabled = false;
            obj.coefsDefault      = zeros(obj.nValidActuator,1);
            obj.coefs             = zeros(obj.nValidActuator,1);
            if isa(obj.modes,'influenceFunction') && ~isempty(p.Results.resolution)
                setInfluenceFunction(obj.modes,obj.nActuator,p.Results.resolution,obj.validActuator);
            end
            obj.log = logBook.checkIn(obj);
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
            display(obj.modes)

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
            if obj.surfaceListener.Enabled
                imagesc(obj);
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