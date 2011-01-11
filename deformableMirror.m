classdef deformableMirror < handle
    % Create a deformableMirror object
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
        coefsUnit = 1;
        % deformableMirror tag
        tag = 'DEFORMABLE MIRROR';
%% Added by Peter Hampton
        % surfaceUnit adjusts the output DM phase to the appropriate
        % units. coefsUnit is an input gain. These are not the same thing.
        % I recommend removing coefsUnit as the implementation can be
        % confusing. If coefsUnit = 1/1000, then dm.coefs = dm.coefs
        % reduces by a factor of 1000 everytime rather than remain static.
        % - PJH
        surfaceUnit = 10^-6;
        
        
        %Distortion Variables
        %system sample time
        sampleTime;
        % distortion based variables 
        % first order time constant ie. 1 - exp(-t/distortionTau)
        distortionTau;          
        % a constant magnitude of distortion
        distortionConstant;     
        % the current state of a filter driving the distortion magnitude
        distortionState;      
        % the seed for the random variable that governs distortion shape
        distortionSeed;         
        % the max dynamic distortion (and is added to constant distortion)
        distortionSaturation;   
        
        %Misalignment Variables
        % the tilt of the DM in the horizontal and vertical direction
        % one unit of tilt has a max phase of +/- 1.9778 microns at the edge when the resolution is 90.
        % Since there are no physical units for mirror width, it is up to
        % the user to determine appropriate magnitude of tilt inputs.
        horizontalTilt;             % magnitude multiplied by normalized tilt mode
        verticalTilt;               % magnitude multiplied by normalized tilt mode
        % the displacement of the DM in the horizontal and vertical directions
        horizontalDisplacement;     % in pixels (can be non-integer)
        verticalDisplacement;       % in pixels (can be non-integer)
        % the rotation of the DM normal to the DM
        rotationAngle;              % in radians
        zern                        % tilt zernike modes

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
        % a flag that enables the distortion code
        distortionOn;
        % the actuation commands that represent distortion
        distortionCoefs;
        % an automatically calculated factor to normalize the 
        % distortion shape to 1 unit so that distortionConstant
        % and distortionSaturation are in the same units as the 
        % phase screen
        distortionFactor;
    end
    
    properties (Access=private)
        p_coefs; 
        p_surface;
        p_validActuator;
        p_distortionOn;
        p_distortionCoefs;
        p_distortionFactor;
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
%-------- pjh - added these parameters in order to deal with actuator drift
            p.addParamValue('distortionSeed',[], @isnumeric); 
            p.addParamValue('distortionTau', inf, @isnumeric);                                    
            p.addParamValue('distortionSaturation', 0, @isnumeric);                      
            p.addParamValue('distortionConstant', 0, @isnumeric);                       
            p.addParamValue('distortionState', 0, @isnumeric); 
            p.addParamValue('sampleTime', 0, @isnumeric);
 %------- pjh - added these parameters to deal with misalignment
            p.addParamValue('horizontalTilt',0, @isnumeric); 
            p.addParamValue('verticalTilt',0, @isnumeric); 
            p.addParamValue('horizontalDisplacement',0, @isnumeric); 
            p.addParamValue('verticalDisplacement',0, @isnumeric);
            p.addParamValue('rotationAngle',0, @isnumeric);
 % ----------------------
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
%-------- pjh - Distortion variables        
            obj.distortionOn        = false;
            obj.distortionSeed      = p.Results.distortionSeed;
            obj.distortionTau       = p.Results.distortionTau;                     
            obj.distortionSaturation= p.Results.distortionSaturation;
            obj.distortionConstant  = p.Results.distortionConstant;
            obj.distortionState     = p.Results.distortionState;       
            obj.sampleTime          = p.Results.sampleTime;
%-----------Misalignment variables
            obj.horizontalTilt      = p.Results.horizontalTilt;
            obj.verticalTilt        = p.Results.verticalTilt;                     
            obj.horizontalDisplacement = p.Results.horizontalDisplacement;
            obj.verticalDisplacement= p.Results.verticalDisplacement;
            obj.rotationAngle       = p.Results.rotationAngle;  
            obj.zern = zernike(1:3,'resolution',p.Results.resolution);
%---------------------
            obj.log = logBook.checkIn(obj);

            display(obj)
        end
        
        %% Destructor
        function delete(obj)
%             if isa(obj.modes,'influenceFunction')
%                 delete(obj.modes)
%             end
            delete(obj.zern)
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
            if obj.distortionOn
                % new distortion based on coefs of last sample period.
                a = exp(-obj.sampleTime./obj.distortionTau); 
                obj.distortionState = (1-a)*obj.distortionSaturation*...
                    (sqrt(obj.coefs(:)'*obj.coefs(:)./obj.nValidActuator)) + a*obj.distortionState;
            end
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
        %% Set and Get distortion coefs
        function out = get.distortionCoefs(obj)
            % if this is the first get then calculate the distortion shape
            % Note that this can be reset by setting the coefs to a scalar
            if isempty(obj.p_distortionCoefs) || isscalar(obj.p_distortionCoefs)
                if ~isempty(obj.distortionSeed)
                    % sets the seed so the distortion is repeatable
                    randn('state',obj.distortionSeed)
                end
                % the gradient of the distortion is random
                dzdx        = randn(obj.nActuator-1);
                dzdy        = randn(obj.nActuator-1);
                [D X Y]     = GradientAnalysis(dzdx,dzdy);
                R           = GradientSynthesisPoisson(D,X,Y,true,true);
                R           = R(1:obj.nActuator,1:obj.nActuator);
                Rmasked     = R.*obj.validActuator;
                Rmean       = sum(sum(Rmasked))/sum(sum(obj.validActuator));
                Rrms        = sqrt(sum(sum((Rmasked-Rmean).^2))./sum(sum(obj.validActuator)));
                R           = (Rmasked-Rmean)./Rrms;
                obj.p_distortionCoefs = R(obj.validActuator == 1);
            end
            out = obj.p_distortionCoefs;
        end
        function set.distortionCoefs(obj,val)
            % allows user to set their own distortion or to set as scalar
            % to reset get routine above.
            obj.p_distortionCoefs = val;
        end
        %% Get and Set the distortion Saturation in nm and gain to convert nm gain to actuation command
        function out = get.distortionFactor(obj)
            % This automatic function approximates the rms with the coarse 
            % pupil defined by valid actuators rather than the high resolution 
            % telescope pupil. The distortionFactor is the inverse of the
            % dm surface rms for the given unadjusted distortionCoefs.
            if isempty(obj.p_distortionFactor)
                obj.p_distortionFactor = 0;
            end
            if obj.p_distortionFactor == 0
                obj.p_distortionFactor = 1; % this is to block recursive loops.
                resolution = size(obj.surface);
                pitch = resolution./(obj.nActuator-1);
                tempSurface = obj.modes*obj.distortionCoefs(:);
                tempSurface = reshape(tempSurface,resolution(1),resolution(2));
                surfaceBin  = zeros(obj.nActuator);
                pixelBin    = zeros(obj.nActuator);
                for row = 0:obj.nActuator-1
                    for col = 0:obj.nActuator-1
                        rS = max([1 floor(1+row*pitch(1)-pitch(1)/2)]);
                        rF = min([resolution(1) floor((1+row)*pitch(1)-pitch(1)/2)]);
                        cS = max([1 floor(1+col*pitch(2)-pitch(2)/2)]);
                        cF = min([resolution(2) floor((1+col)*pitch(2)-pitch(2)/2)]);
                        surfaceBin(1+row,1+col) = sum(sum(tempSurface(rS:rF,cS:cF).*tempSurface(rS:rF,cS:cF)));
                        pixelBin(1+row,1+col)   = (rF-rS+1)*(cF-cS+1);
                    end
                end
                maskedSurface = surfaceBin.*obj.validActuator;
                rmsSurface = sqrt(sum(maskedSurface(:))/sum(pixelBin(:).*obj.validActuator(:)));
                obj.p_distortionFactor = 1/rmsSurface;
            end
            out = obj.p_distortionFactor;
        end
        function set.distortionFactor(obj,val)
            % If the user calculates a better distortion Factor that
            % results in an rms of 1 when distortionSaturation = 1, and the
            % rms of the coefs = 1 then they can use this. The automatic
            % function approximates the rms with the coarse pupil defined
            % by valid actuators rather than the high resolution telescope
            % pupil.
            obj.p_distortionFactor = val;
            fprintf('Distortion factor will be recalculated when set to 0\n')
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
            % modified this function from a simple coeficient projection to
            % include actuator drift due to current in the coils.
            
            % The distortion is defined to have a saturation of
            % obj.distortionGain in the DM surface units when the obj.coefs is 1 rms
            % The obj.distortionConstant is a constant magnitude
            % of the distortion. The magnitude of the distortion is always
            % positive by this definition (unless Gain or Constant is defined as negative).
            
            % The sample time should match the simulation so that a simple
            % definition of tau will define the reaction time of the
            % distortion. tau is the time constant of a first order system.
            % ie. 1 - exp(-t/tau) defines a first order unit step response.
            
            % A major assumption that may be sufficient for simulation (only
            % if it is not corrected) is that the distortion pattern does
            % not change depending on surface shape and that the magnitude
            % is only dependant on the rms of the coefs (which loosely
            % models the rms of the current in the DM coils. Improvements 
            % of this can be accommodated by changing the distortionCoefs.
            
            % distortionOn was included to enable/disable this code with
            % out worrying about whether all the coefficients are properly
            % set/cleared. Leave distortionOn = false for interaction
            % matrices because the time scale is in minutes and this function 
            % only works for coomand vectors, not command matrices.
            
            if obj.distortionOn
                obj.p_surface = obj.modes*(obj.coefs(:)...
                    + obj.distortionCoefs(:)*obj.distortionFactor*(obj.distortionState + obj.distortionConstant));
            else
                obj.p_surface = obj.modes*obj.coefs; % this line was out of this if statement originally
            end
          
            if obj.lex
                dmShape = obj.p_surface;
            else
                dmShape = utilities.toggleFrame(obj.p_surface,3);
            end
        end
                %% Observe On state
        function out = get.distortionOn(obj)
             if isempty(obj.p_distortionOn)
                 % assume Distortion is OFF if empty
                 obj.p_distortionOn = false;
             end
             out = obj.p_distortionOn;
             
        end
        %% Turn Distortions On/Off
        function set.distortionOn(obj,val)
            if val
                dot_dC = obj.distortionCoefs'*obj.distortionCoefs;
                if dot_dC > 0 && obj.distortionTau > 0 && obj.distortionTau < inf && obj.distortionSaturation > 0
                    obj.p_distortionOn = true;
                    fprintf('Mirror will distort on each sample.\n')
                    return
                else
                    obj.p_distortionOn = false;
                    fprintf('The distortion is not fully defined! Mirror will not distort!\n')
                end
            else
                obj.p_distortionOn = false;
                obj.distortionState = 0;
                fprintf('The distortion is disabled and reset.\n')
            end
        end
        function relay(obj,src)
            %% RELAY deformable mirror to source relay
            %
            % relay(obj,srcs) writes the deformableMirror amplitude and
            % phase into the properties of the source object(s)
            
            nSrc       = numel(src);
            wavenumber = 2*pi/src(1).wavelength;
            dmPhase    = obj.surface;
            % Misalignments ------------------
            if obj.horizontalTilt ~= 0 || obj.verticalTilt ~= 0
                dmPhase = dmPhase +...
                        obj.horizontalTilt*reshape(obj.zern.modes(:,2),size(dmPhase)) +...
                        obj.verticalTilt*reshape(obj.zern.modes(:,3),size(dmPhase));
            end
            dmPhase = -2*dmPhase*obj.surfaceUnit*wavenumber; %original code was dmPhase = -2*obj.surface*wavenumber
            if obj.rotationAngle ~= 0 || obj.horizontalDisplacement ~= 0 || obj.verticalDisplacement ~= 0
                dmPhase = rotateDisplace(dmPhase, obj.rotationAngle, obj.horizontalDisplacement, obj.verticalDisplacement);
                if size(src.amplitude(:),1) == 1                    
                    src.amplitude = src.amplitude*ones(size(dmPhase));
                end
                src.amplitude(isnan(dmPhase)) = 0;
                dmPhase(isnan(dmPhase)) = 0;
            end


            % Currently assuming perfect alignment in its position along
            % optical axis. zLocation could probably be used for optical
            % axis misalignment.
            % --------------
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
            
            if isa(obj.modes,'zernike')
                out = zernikeStats.residualVariance(obj.modes.nMode,telAtm.opticalAberration,telAtm);
            else
                d = telAtm.D/(obj.nActuator-1);
                fc = 1/d/2;
                out = phaseStats.variance(telAtm.opticalAberration) - ...
                    dblquad( @(fx,fy) phaseStats.spectrum( hypot(fx,fy) , telAtm.opticalAberration ) , ...
                    -fc,fc,-fc,fc);
            end
            if nargin>2
                if nargin<4
                    unit = 1; % Shouldn't this be 0 to default to meters? - PJH
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