classdef deformableMirror < handle
   
    properties
        % # of actuator
        nActuator;
        % influence functions
        modes;
        % electronic driver
        driver;
        % optical aberrations
        opticalAberration;
        % coefficients default
        coefsDefault;
        % valid actuator
        validActuator;
        % surface listener
        surfaceListener;
    end
    
    properties (SetObservable=true,Dependent,SetAccess=private)
        % the shape of the DM
        surface;
    end
    
    
    properties (Dependent)
        % # of actuators in the pupil
        nValidActuator;
        % coefficients
        coefs;
    end
    
    properties (Access=private)
        p_coefs; 
        imageHandle;
    end
        
    methods
        
        % Constructor
        function obj = deformableMirror(nActuator,varargin)
            p = inputParser;
            p.addRequired('nActuator', @isnumeric);
            p.addParamValue('modes', [], @(x) isnumeric(x) || isa(x,'influenceFunction') );
            p.addParamValue('resolution', [], @isnumeric);
            p.addParamValue('validActuator', ones(nActuator), @islogical);
            p.addParamValue('opticalAberration', [], @(x) isnumeric(x) || isa(x,'telescope') || isa(x,'zernike'));
            p.parse(nActuator, varargin{:});
            obj.nActuator         = p.Results.nActuator;
            obj.validActuator     = p.Results.validActuator;
            obj.modes             = p.Results.modes;
            obj.surfaceListener = addlistener(obj,'surface','PostSet',...
                @(src,evnt) obj.imagesc );
            obj.surfaceListener.Enabled = false;
            obj.coefsDefault      = zeros(obj.nValidActuator,1);
            obj.coefs             = zeros(obj.nValidActuator,1);
            obj.opticalAberration = p.Results.opticalAberration;
            if isa(obj.modes,'influenceFunction')
                setInfluenceFunction(obj.modes,obj.nActuator,p.Results.resolution,obj.validActuator);
            end
        end
        
%         % Destructor
%         function delete(obj)
%             if isa(obj.modes,'influenceFunction')
%                 delete(obj.modes)
%             end
%         end
        
        function out = get.nValidActuator(obj)
            out = sum(obj.validActuator(:));
        end
        
        % Set and Get coefs
        function out = get.coefs(obj)
            out = obj.p_coefs;
        end
        function set.coefs(obj,val)
            obj.p_coefs = val + obj.coefsDefault;
            if isa(obj.driver,'function_handle')
                obj.driver(obj.p_coefs);
                return
            end
            if obj.surfaceListener.Enabled
                imagesc(obj);
            end
        end
        
        % Get the dm shape
        function dmShape = get.surface(obj)
            dmShape = obj.modes*obj.coefs;
            if isempty(obj.opticalAberration)
                dmShape = utilities.toggleFrame(dmShape,3);
            else
                switch class(obj.opticalAberration)
                    case 'double'
                        dmShape = utilities.toggleFrame(dmShape,3);
                    case {'telescope','zernike'}
                        dmShape(~obj.opticalAberration.pupil,:) = 0;
                        dmShape = utilities.toggleFrame(dmShape,3);
                    otherwise
                end
            end
        end
        
        function dmShape = getPhase(obj,src)
            % GETPHASE Shape of the deformable mirror
            %
            % dmShape = getPhase(obj) computes the DM shape. If
            % opticalAberration property is set, the DM shape is
            % substracted from the opticalAberration
            %
            % dmShape = getPhase(obj,src) passes the source to
            % opticalAberration if needed
            
            if nargin<2
                src = [];
            end
            dmShape = obj.modes*obj.coefs;
            if isempty(obj.opticalAberration)
                dmShape = utilities.toggleFrame(dmShape,3);
            else
                switch class(obj.opticalAberration)
                    case 'double'
                        dmShape = utilities.toggleFrame(dmShape,3);
                        dmShape = obj.opticalAberration - dmShape*2;
                    case 'telescope'
                        dmShape(~obj.opticalAberration.pupil,:) = 0;
                        dmShape = utilities.toggleFrame(dmShape,3);
                        phaseScreen = getPhaseScreen(obj.opticalAberration,src);
                        [n,m] = size(phaseScreen);
                        dmShape = repmat(dmShape,1,m/n);
                        dmShape = bsxfun(@minus,phaseScreen,dmShape*2);
                    case 'zernike'
                        dmShape(~obj.opticalAberration.pupil,:) = 0;
                        dmShape = utilities.toggleFrame(dmShape,3);
                        dmShape = bsxfun(@minus,obj.opticalAberration.phase,dmShape*2);
                    otherwise
                end
            end
        end
        
        function dmShape = getPistonFreePhase(obj,src)
            % GETPISTONFREEPHASE Shape of the deformable mirror
            %
            % dmShape = getPhase(obj) computes the DM shape. If
            % opticalAberration property is set, the DM shape is
            % substracted from the opticalAberration
            %
            % dmShape = getPhase(obj,src) passes the source to
            % opticalAberration if needed
            %
            % The piston mode is removed from dmShape
            %
            % See also deformableMirror/getPhase
            
            if nargin<2
                src = [];
            end
            dmShape = getPhase(obj,src);
            u = obj.opticalAberration.pupilLogical;
            dmShape(u) = dmShape(u) - mean(dmShape(u));
            
        end
        
        function out = getWave(obj,src)
            if nargin<2
                src = [];
            end
            out = exp(1i.*getPhase(obj,src));
            if isa(obj.opticalAberration,'telescope') || isa(obj.opticalAberration,'zernike')
                out = utilities.toggleFrame(out,2);
                out(~obj.opticalAberration.pupil,:) = 0;
                out = utilities.toggleFrame(out,3);
            end
        end
        
        function varargout = imagesc(obj,varargin)
            % IMAGESC Display the lenslet imagelets
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