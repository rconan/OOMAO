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
            p.addParamValue('modes', [], @(x) isnumeric(x) || ...
                (isa(x,'influenceFunction') || isa(x,'zernike') || isa(x,'hexagonalPistonTipTilt')) );
            p.addParamValue('resolution', [], @isnumeric);
            p.addParamValue('validActuator', ones(nActuator), @islogical);
            p.addParamValue('zLocation', 0, @isnumeric);
            p.addParamValue('offset', [0,0], @isnumeric);
            p.parse(nActuator, varargin{:});
            obj.nActuator         = p.Results.nActuator;
            obj.p_validActuator     = p.Results.validActuator;
            obj.modes             = p.Results.modes;
            obj.zLocation             = p.Results.zLocation;
            setSurfaceListener(obj)
            if ( isa(obj.modes,'influenceFunction') || isa(obj.modes,'hexagonalPistonTipTilt')) && ~isempty(p.Results.resolution)
                setInfluenceFunction(obj.modes,obj.nActuator,...
                    p.Results.resolution,obj.validActuator,1,p.Results.offset);
            elseif isa(obj.modes,'zernike')
                obj.p_validActuator = true(1,obj.modes.nMode);
            end
            if isa(obj.modes,'hexagonalPistonTipTilt')
                obj.coefsDefault      = zeros(3*obj.nValidActuator,1);
                obj.coefs             = zeros(3*obj.nValidActuator,1);
            else
                obj.coefsDefault      = zeros(obj.nValidActuator,1);
                obj.coefs             = zeros(obj.nValidActuator,1);
            end
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
        
        function obj = saveobj(obj)
            %% SAVEOBJ
            delete(obj.surfaceListener)
            add(obj.log,obj,'Save!')
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
        
        
        function obj = mldivide(obj,src)
            %% \ Least square fit to the influence functions
            %
            % obj = obj\src fits the source object wavefront onto the
            % deformable mirror object influence functions and stro the
            % projection coefficients into the object coefficients vector

            F = obj.modes.modes(src.mask,:);
            maps = utilities.toggleFrame(src.phase,2);
            obj.coefs = 0.5*(F\maps(src.mask,:))/src.waveNumber;
        end
        
        function calib = calibration(obj,sensor,src,calibDmStroke,steps)
            %% CALIBRATION DM calibration
            %
            % obj = calibration(obj,sensor,src,calibDmCommands) calibrate
            % the DM object with the sensor object using the calibration
            % source src and the actuator stroke calibDmStroke
            %
            % obj = calibration(obj,sensor,src,calibDmCommands,nSteps)
            % calibrate the DM object with the sensor object using the
            % calibration source src and the actuator stroke calibDmStroke;
            % the calibration process is split in nSteps steps.
            
            obj.coefs = 0;
            if nargin<5
                steps = 1;
                if obj.nValidActuator>1000
                    steps  = 25;
                end
            end
            
            if sensor.lenslets.nLenslet>1
                
                src = src*obj*sensor;
                
                if isscalar(calibDmStroke)
                    calibDmCommands = speye(obj.nValidActuator)*calibDmStroke;
                else
                    calibDmCommands = calibDmStroke;
                    calibDmStroke = 1;
                end
                
                tId = tic;
                if steps==1
                    obj.coefs = calibDmCommands;
                    +src;
                    pokeMatrix = sensor.slopes;
                else
                    nMode = size(calibDmCommands,2);
                    nC              = floor(nMode/steps);
                    u               = 0;
                    pokeMatrix  = zeros(sensor.nSlope,nMode);
%                     fprintf(' . actuators range:          ')
                    h = waitbar(0,'DM/WFS calibration ...');
                    while u(end)<nMode
                        u = u(end)+1:min(u(end)+nC,nMode);
                        waitbar(u(end)/nMode)
%                         fprintf('\b\b\b\b\b\b\b\b\b%4d:%4d',u(1),u(end))
                        obj.coefs = calibDmCommands(:,u);
                        +src;
                        pokeMatrix(:,u) = sensor.slopes;
                    end
                    close(h)
                end
                elt = toc(tId);
                
%                 pokeMatrix = src.wavelength*pokeMatrix./calibDmStroke;
                fprintf('__ Poke Matrix Stats ___\n')
                fprintf(' . computing time: %5.2fs\n',elt)
                fprintf(' . size: %dx%d\n',size(pokeMatrix))
                nonZerosIndex = pokeMatrix(:)~=0;
                nonZeros = sum(nonZerosIndex);
                fprintf(' . non zeros values: %d i.e. %4.2f%%\n',nonZeros,100*nonZeros/numel(pokeMatrix))
                pokeMatrixNZ = pokeMatrix(nonZerosIndex);
                fprintf(' . min. and max. values: [%5.2f,%5.2f]\n',max(pokeMatrixNZ),min(pokeMatrixNZ))
                fprintf(' . mean and median of absolute values: [%5.2f,%5.2f]\n',mean(abs(pokeMatrixNZ)),median(abs(pokeMatrixNZ)))
                fprintf('________________________\n')
                pokeMatrix = pokeMatrix./calibDmStroke;
                calib = calibrationVault(pokeMatrix);
                
            else 
                
                %%% Tip-Tilt sensor calibration
                
                tel = src.opticalPath{1};
                zern = zernike(tel,2:3);
                zern.c = eye(2)*src.wavelength/4;
                src = src.*zern;
                buf = reshape(src.phase,tel.resolution,2*tel.resolution);
                
                % zernike projection onto DM influence functions
                obj = obj\src;
                src = src.*tel*obj;
                dmTtCoefs = obj.coefs;
                
                buf = [buf;reshape(src.phase,tel.resolution,2*tel.resolution)];
                figure
                imagesc(buf)
                axis square
                colorbar
                
                obj.coefs = dmTtCoefs*calibDmStroke;
                src = src.*tel*obj*sensor;
                pokeTipTilt = sensor.slopes/calibDmStroke;
                calib = calibrationVault(pokeTipTilt);
                calib.spaceJump = dmTtCoefs;
                calib.nThresholded = 0;
            end
            
            obj.coefs = 0;
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
            switch class(obj.modes)
                case 'zernike'
                    out = zernikeStats.residualVariance(obj.modes.nMode,atm,telAtm);
                case 'influenceFunction'
                    d = telAtm.D/(obj.nActuator-1);
                    fc = 1/d/2;
                    a = phaseStats.variance(atm);
                    b = dblquad( @(fx,fy) phaseStats.spectrum( hypot(fx,fy) , atm ) , ...
                        -fc,fc,-fc,fc);
                    out = a - b;
                case 'hexagonalPistonTipTilt'
                    nCycle = roots([3,3,1-obj.nActuator/3]);
                    nCycle(nCycle<0) = [];
                    smallTel = telescope(telAtm.D*0.5/(nCycle-1));
                    out = zernikeStats.residualVariance(3,atm,smallTel);                    
            end
            atm.wavelength = atmWavelength;
            if nargin>3
                out = 10^-unit*sqrt(out)/src.waveNumber;
            end            
        end
        
        function varargout = imagesc(obj,varargin)
            %% IMAGESC Display the deformable mirror surface
            %
            % imagesc(obj) displays the deformable mirror surface
            %
            % imagesc(obj,'PropertyName',PropertyValue) displays the
            % deformable mirror surface and set the properties of the
            % graphics object imagesc
            %
            % h = imagesc(obj,...) returns the graphics handle
            %
            % See also: imagesc
            
            if ishandle(obj.imageHandle)
                set(obj.imageHandle,'Cdata',obj.surface*1e6,varargin{:});
            else
                %                 figure
                obj.imageHandle = image(obj.surface*1e6,...
                    'CDataMApping','Scaled',varargin{:});
                set( get(obj.imageHandle,'parent') ,'visible','off')
                colormap(pink)
                axis square
                ylabel(colorbar,'[micron]')
            end
            if nargout>0
                varargout{1} = obj.imageHandle;
            end
        end
        
        
    end
        
    methods (Static)
            
        function obj = loadobj(obj)
            %% LOADOBJ
            add(obj.log,obj,'Load!')
            setSurfaceListener(obj)
            obj.log = logBook.checkIn(obj);
        end
        
    end
    
    methods (Access=private)
        
        function setSurfaceListener(obj)
            %% SETSLOPESLISTENER Slopes listener
            obj.surfaceListener = addlistener(obj,'surface','PostSet',...
                @(src,evnt) obj.imagesc );
            obj.surfaceListener.Enabled = false;
        end
        
    end

end