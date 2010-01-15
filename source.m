classdef source < stochasticWave
    % SOURCE Create a source object
    %
    % src = source creates an on-axis star at infinity
    %
    % src = source('parameter',value) creates a source object from
    % parameter-value pair arguments. The parameters are zenith, azimuth,
    % height, wavelength, magnitude, nPhoton, width and asterism.
    %
    % Example:
    % src = source;
    % creates an on-axis source
    % src =
    % source('zenith',30*cougarConstants.arcsec2radian,'azimuth',pi/4); 
    % creates an off-axis source at 30 arcsec zenith and 45degree azimuth
    % src = source('asterism',{[0,0],[5,60*cougarConstants.arcsec2radian,pi/3]})
    % set an asterism with one star on-axis and five star on a 60 arcsec
    % radius shifted of 30 degrees
    
    properties
        % source zenith angle
        zenith;
        % source azimuth angle
        azimuth;
        % source height
        height;
        % source wavelength
        wavelength;
        % source full-width-half-max
        width;
        % source view point
        viewPoint;
        % # of photon [m^{-2}.s^{-1}] 
        nPhoton;
        % cell array of handles of objects the source is propagating through  
        opticalPath;
    end
    
    properties (SetAccess=private)
        directionVector;
        wavefront;
    end
    
    properties (Dependent)
        waveNumber;
        magnitude;
    end
    
    properties (Access=private)
%         p_nPhoton;
        p_magnitude;
        tel;
        log;
    end
        
    methods
        
        %% Constructor
        function obj = source(varargin)
            obj = obj@stochasticWave;
            p = inputParser;
            p.addParamValue('asterism',[],@iscell);
            p.addParamValue('zenith',0,@isnumeric);
            p.addParamValue('azimuth',0,@isnumeric);
            p.addParamValue('height',Inf,@isnumeric);
            p.addParamValue('wavelength',[],@isnumeric);
            p.addParamValue('magnitude',[],@isnumeric);
            p.addParamValue('nPhoton',[],@isnumeric);
            p.addParamValue('width',0,@isnumeric);
            p.addParamValue('viewPoint',[0,0],@isnumeric);
            p.parse(varargin{:});
            nAst    = numel(p.Results.asterism);
            nHeight = numel(p.Results.height);
            nObj = nAst + nHeight;
            if nObj>1
                if nAst>0
                    ast = p.Results.asterism;
                else
                    ast = {[p.Results.zenith,p.Results.azimuth]};
                nAst = nAst + 1;
                end
                z = zeros(1,nAst);
                a = zeros(1,nAst);
                nObj = 0;
                for kAst = 1:nAst
                    if length(ast{kAst})==3
                        n = ast{kAst}(1);
                        z(nObj+1:nObj+n) = ast{kAst}(2).*ones(1,n);
                        a(nObj+1:nObj+n) = ast{kAst}(3) + (0:n-1)*2*pi/n;
                        nObj = nObj + n;
                    else
                        z(nObj+1) = ast{kAst}(1);
                        a(nObj+1) = ast{kAst}(2);
                        nObj = nObj + 1;
                    end
                end
                obj( 1 , nObj , nHeight ) = source;
                for kObj = 1:nObj
                    for kHeight = 1:nHeight
                        obj(1,kObj,kHeight).zenith     = z(kObj);
                        obj(1,kObj,kHeight).azimuth    = a(kObj);
                        obj(1,kObj,kHeight).height     = p.Results.height(kHeight);
                        obj(1,kObj,kHeight).wavelength = p.Results.wavelength;
                        obj(1,kObj,kHeight).nPhoton    = p.Results.nPhoton;
                        obj(1,kObj,kHeight).magnitude  = p.Results.magnitude;
                        obj(1,kObj,kHeight).width      = p.Results.width;
                        obj(1,kObj,kHeight).viewPoint  = p.Results.viewPoint;
                        setDirectionVector(obj(1,kObj,kHeight))
                    end
                end
            else
                obj.zenith     = p.Results.zenith;
                obj.azimuth    = p.Results.azimuth;
                obj.height     = p.Results.height;
                obj.wavelength = p.Results.wavelength;
                obj.nPhoton    = p.Results.nPhoton;
                obj.magnitude  = p.Results.magnitude;
                obj.width      = p.Results.width;
                obj.viewPoint  = p.Results.viewPoint;
                setDirectionVector(obj)
            end
            obj.log = logBook.checkIn(obj);
        end
        
        %% Destructor
        function delete(obj)
            checkOut(obj.log,obj)
        end
        
        %% Get and Set nPhoton
%         function set.nPhoton(obj,val)
%             obj.p_nPhoton = val;
%         end
        function out = get.nPhoton(obj)
            if isempty(obj.nPhoton)
                if ~isempty(obj.p_magnitude)
                    index = abs(photometry.wavelengths-photometry.R)<=eps(max(photometry.wavelengths));
                    out = photometry.deltaWavelengths(index).*...
                        1e6.*photometry.wavelengths(index).*photometry.e0(index).*...
                        10.^(-0.4*obj.p_magnitude)./(cougarConstants.plank*cougarConstants.c);
                else
                    out = [];
                end
            else
                out = obj.nPhoton;
            end
        end
        %% Get and Set magnitude
        function out = get.magnitude(obj)
            out = obj.p_magnitude;
        end
        function set.magnitude(obj,val)
            obj.nPhoton = [];
            obj.p_magnitude = val;
        end

        
%         function bool = eq(obj1,obj2)
%             % == (EQ) Sources comparison
%             % src1==src2 returns true if both objects have the same zenith
%             % and azimuth angles and the same height
%             
%             bool = false;
%             if obj1.zenith==obj2.zenith && ...
%                     obj1.azimuth==obj2.azimuth && ...
%                     obj1.height==obj2.height
%                 bool = true;
%             end
%         end
        
        function bool = isNgs(obj)
            %% ISNGS NGS check
            % bool = isNgs(obj) returns true if the height of the source
            % object is infinite
            
            bool = isinf(obj.height);
        end
        
        function bool = isLgs(obj)
            %% ISLGS LGS check
            % bool = isLgs(obj) returns true if the height of the source
            % object is finite
            
            bool = isfinite(obj.height);
        end
        
        function bool = isOnAxis(obj)
            %% ISONAXIS On-axis check
            % bool = isOnAxis(obj) returns true if both zenith and azimuth
            % angles are zero
            
            bool = (obj.zenith==0) && (obj.azimuth==0);
        end
        
        function setDirectionVector(obj)
            obj.directionVector.x = tan(obj.zenith)*cos(obj.azimuth);
            obj.directionVector.y = tan(obj.zenith)*sin(obj.azimuth);
            obj.directionVector.z = 1;
        end
        
        function val = get.waveNumber(obj)
            val = 2*pi/obj.wavelength;
        end        
        
        function varargout = reset(obj)
            %% RESET Reset wave properties
            %
            % reset(obj) resets the mask to [], the amplitude to 1, the
            % phase to 0 and the optical path to [];
            %
            % obj = reset(obj) resets and returns the reset object
            
            for kObj = 1:numel(obj);
                obj(kObj).mask        = [];
                obj(kObj).p_amplitude = 1;
                obj(kObj).p_phase     = 0;
                obj(kObj).opticalPath = [];
            end
            if nargout>0
                varargout{1} = obj;
            end
        end

        function obj = ge(obj,otherObj)
            %% >= Source object propagation operator
            %
            % src = src>=otherObj propagate src through otherObj
            % multiplying the source amplitude by the otherObj transmitance
            % and adding the otherObj phase to the source phase
            
            if isa(otherObj,'shackHartmann')
                propagateThrough(otherObj.lenslets,obj)
                grabAndProcess(otherObj)
            elseif isa(otherObj,'lensletArray')
                propagateThrough(otherObj,obj)
            elseif ~isempty(otherObj) % do nothing if the other object is empty
                nObj = numel(obj);
                if nObj>1 % if the source is an array, recursive self-calls
                    for kObj = 1:nObj
                        ge(obj(kObj),otherObj);
                    end
                else
                    if isobject(otherObj)
                        relay(otherObj,obj);
                    else
                        obj.amplitude = abs(otherObj);
                        obj.phase     = angle(otherObj);
                    end
                end
            end
            obj.opticalPath{ length(obj.opticalPath)+1 } = otherObj;
        end
        
        function obj = eq(obj,otherObj)
            %% == Source object propagation operator
            %
            % src==otherObj propagate src through otherObj setting the
            % source amplitude to the otherObj transmitance and the source
            % phase to the otherObj phase
            %
            % src = src==otherObj returns the source object
            
             ge(reset(obj),otherObj);
        end

        function out = fresnelPropagation(obj,tel)
            %% FRESNELPROPAGATION Source propagation to the light collector
            % fresnelPropagation(a,tel) propagates the source seen
            % from the given view point to the telescope primary mirror
            % out fresnelPropagation(a,tel) propagates the source seen
            % from the given view point to the telescope primary mirror and
            % returns the wavefront in radian

            if isempty(obj.tel) || ~isequal(obj.tel,tel)
                fprintf(' @(source)> Computing the wavefront ...\n')
                obj.tel = tel;
                if ( numel(obj.height)==1 && isinf(obj.height) ) || ...
                        ( numel(obj.height)==1 && obj.height==obj.tel.focalDistance )                  
                    obj.wavefront = zeros(tel.resolution);
                else
                    rho     = utilities.cartAndPol(tel.resolution,tel.R,...
                        'offset',obj.viewPoint,'output','radius');
                    if isinf(tel.focalDistance)
                        s0 = 0;
                    else
                        s0 = hypot(rho,tel.focalDistance);
                    end
                    h       = reshape(obj.height,[1,1,length(obj.height)]);
                    s = sqrt(bsxfun(@plus,rho.^2,h.^2));
                    obj.wavefront = bsxfun(@minus,s,s0);
                    obj.wavefront = 2*pi*obj.wavefront/obj.wavelength;
                    % 2\pi demodulation for a meamingfull phase
                    obj.wavefront = mod(obj.wavefront,2*pi);
                    %                     obj.wavefront = exp(1i.*2*pi*obj.wavefront/obj.wavelength)./s;
                end
            end
            if nargout>0
                out = obj.wavefront;
            end
        end
        
        function copy = clone(obj)
            %% CLONE Create an object clone
            %
            % copy = clone(obj) clones the source object to another object
            % with a different handle
            
            copy = eval(class(obj));
            meta = eval(['?',class(obj)]);
            for p = 1: size(meta.Properties,1)
                pname = meta.Properties{p}.Name;
                try
                    if ~meta.Properties{p}.Dependent
                        eval(['copy.',pname,' = obj.',pname,';']);
                    end
                catch ME
                    fprintf(['\nCould not copy ',pname,'.\n']);
                    rethrow(ME)
                end
            end
        end
        
    end
    
end