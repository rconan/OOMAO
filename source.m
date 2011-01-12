classdef source < stochasticWave & hgsetget
    %% SOURCE Create a source object
    %
    % src = source creates an on-axis star at infinity
    %
    % src = source('parameter',value) creates a source object from
    % parameter-value pair arguments. The parameters are zenith, azimuth,
    % height, wavelength, magnitude, nPhoton, width and asterism.
    %
    % The source class represents celestial objects. It inherits from
    % stochasticWave which contains the phase and amplitude of the wave
    % emitted by the source object. The source propagates from one object
    % to another object on the optical path by using the propagation
    % operator ".*" and "*". The propagation through a pre-set optical
    % path can be re-played with the unary plus (+obj) operator.
    %
    % Example :
    % To create an on-axis source: src = source; 
    %
    % To create an off-axis source at 30 arcsec zenith and 45degree
    % azimuth: 
    % src = source('zenith',30*constants.arcsec2radian,'azimuth',pi/4); 
    %
    % To set an asterism with one star on-axis and five star on a 60 arcsec
    % radius shifted of 30 degrees : 
    % src = source('asterism',{[0,0],[5,60*constants.arcsec2radian,pi/3]})
    %
    % Setting the magnitude
    % src = source('zenith',30*constants.arcsec2radian,'wavelength',photometry.H,'magnitude',12)
    %
    % For an asterism, all the sources have the same wavelength but they
    % may have different magnitudes
    % src = source('asterism',{[0,0],[5,60*constants.arcsec2radian,0]},'wavelength',photometry.H,'magnitude',[8 10 12 9 11 14])
    %
    % See also: constants, photometry, telescope and atmosphere
    
    properties
        % source height
        height;
        % source full-width-half-max
        width;
        % # of photon [m^{-2}.s^{-1}] 
        nPhoton;
        % cell array of handles of objects the source is propagating through  
        opticalPath;
        log;
        % source tag
        tag;
        % source #
        nSrc;
        % time stamp
        timeStamp = 0;
    end
    
    properties (SetAccess=private)
        directionVector = zeros(3,1);
        wavefront;
    end
    
    properties (Dependent)
        % source zenith angle
        zenith;
        % source azimuth angle
        azimuth;
        % source wavenumber
        waveNumber;
        % source magnitudde
        magnitude;
        % source wavelength
        wavelength;
        % source view point
        viewPoint;
    end
    
    properties (Dependent,SetAccess=private)
        wavelengthInMicron;
    end
    
    properties (Access=private)
%         p_nPhoton;
        p_zenith;
        p_azimuth;
        p_magnitude;
        photometry;
        p_viewPoint;
        tel;
    end
        
    methods
        
        function obj = source(varargin)

            obj = obj@stochasticWave;
            p = inputParser;
            p.addParamValue('asterism',[],@iscell);
            p.addParamValue('zenith',0,@isnumeric);
            p.addParamValue('zenithInArcsec',[],@isnumeric);
            p.addParamValue('zenithInArcmin',[],@isnumeric);
            p.addParamValue('azimuth',0,@isnumeric);
            p.addParamValue('height',Inf,@isnumeric);
            p.addParamValue('wavelength',photometry.V,@(x) isa(x,'photometry'));
            p.addParamValue('magnitude',0,@isnumeric); % Vega magnitude (default)
            p.addParamValue('nPhoton',[],@isnumeric);
            p.addParamValue('width',0,@isnumeric);
            p.addParamValue('viewPoint',[0,0],@isnumeric);
            p.addParamValue('tag','SOURCE',@ischar);
            p.parse(varargin{:});
            persistent nCall
%             if ~isempty(p.Results.zenithInArcsec)
%                 obj.zenith      = p.Results.zenithInArcsec./cougarConstants.radian2arcsec;
%             elseif ~isempty(p.Results.zenithInArcmin)
%                 obj.zenith      = p.Results.zenithInArcmin./cougarConstants.radian2arcmin;
%             end
            if nargin>0 || isempty(nCall)
                nCall = 1;
                nAst    = numel(p.Results.asterism);
                if nAst>0
                    ast = p.Results.asterism;
                else
                    n = length(p.Results.zenith);
                    ast = mat2cell([p.Results.zenith(:),p.Results.azimuth(:)],ones(n,1),2);
                    nAst = nAst + n;
                end
                nHeight = numel(p.Results.height);
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
                magnitude = p.Results.magnitude;
                nMag = length(magnitude);
                obj( 1 , nObj , nHeight ) = source;
                for kObj = 1:nObj
                    for kHeight = 1:nHeight
                        obj(1,kObj,kHeight).p_zenith     = z(kObj);
                        obj(1,kObj,kHeight).p_azimuth    = a(kObj);
                        obj(1,kObj,kHeight).height     = p.Results.height(kHeight);
                        obj(1,kObj,kHeight).wavelength = p.Results.wavelength;
                        obj(1,kObj,kHeight).nPhoton    = p.Results.nPhoton;
                        if ~isempty(magnitude)
                            obj(1,kObj,kHeight).magnitude  = p.Results.magnitude(min(nMag,kObj));
                        end
                        obj(1,kObj,kHeight).width      = p.Results.width;
                        obj(1,kObj,kHeight).viewPoint  = p.Results.viewPoint;
                        obj(1,kObj,kHeight).tag        = p.Results.tag;
                        obj(1,kObj,kHeight).nSrc       = nObj*nHeight;
                        setDirectionVector(obj(1,kObj,kHeight))
                   end
                end
                nCall = [];
            else
                obj.p_zenith     = p.Results.zenith;
                obj.p_azimuth    = p.Results.azimuth;
                obj.height     = p.Results.height;
                obj.wavelength = p.Results.wavelength;
                obj.nPhoton    = p.Results.nPhoton;
                obj.magnitude  = p.Results.magnitude;
                obj.width      = p.Results.width;
                obj.viewPoint  = p.Results.viewPoint;
                obj.tag        = p.Results.tag;
                obj.nSrc       = 1;
                setDirectionVector(obj)
            end
            % Here we are registering the source object(s) into the log
            % book
            if isempty(nCall)
                set(obj,'log',logBook.checkIn(obj));
                display(obj)
            end
         end
        
        %% Destructor
        function delete(obj)
            persistent countObj
            if isempty(countObj)
                countObj = obj.nSrc;
            end
            if ~isempty(obj.log)
                countObj = countObj - 1;
                if countObj==0
                    obj.log.verbose = true;
                    countObj = [];
                else
                    obj.log.verbose = false;
                end
                checkOut(obj.log,obj)
            end
        end
        
        %% Get and Set zenith
        function out = get.zenith(obj)
            out = obj.p_zenith;
        end
        function set.zenith(obj,val)
            obj.p_zenith = val;
            setDirectionVector(obj)
        end
        
        %% Get and Set azimuth
        function out = get.azimuth(obj)
            out = obj.p_azimuth;
        end
        function set.azimuth(obj,val)
            obj.p_azimuth = val;
            setDirectionVector(obj)
        end
        
        %% Get and Set magnitude
        function out = get.magnitude(obj)
            out = obj.photometry.magnitude;
        end
        function set.magnitude(obj,val)
            obj.photometry.magnitude = val;
            if ~isempty(obj.photometry)
                obj.nPhoton = obj.photometry.nPhoton;
                fprintf(' @(source)> # of photon m^{-2}.s^{-1}: %4.2f\n',obj.nPhoton)
            end
        end
        
        %% Get the wavelength in micron
        function out = get.wavelengthInMicron(obj)
            out = obj.photometry.wavelength*1e6;
        end
        function out = get.wavelength(obj)
            out = obj.photometry.wavelength;
        end
        function set.wavelength(obj,val)
            if ~isa(val,'photometry')
                error('oomao:source:wavelength','The wavelength must be set with the photometry class!')
            end
            obj.photometry = val;
            if ~isempty(obj.photometry.magnitude)
                obj.nPhoton = obj.photometry.nPhoton;
            end
        end
        
        %% Get the wavelength in micron
        function out = get.viewPoint(obj)
            out = obj.p_viewPoint;
        end
        function set.viewPoint(obj,val)
            obj.p_viewPoint = val;
            obj.tel         = [];
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
        
        function display(obj)
            %% DISP Display object information
            %
            % disp(obj) prints information about the source object
            [n,m,k] = size(obj);
            nObj = n*m;
            fprintf('___ %s ___\n',obj(1).tag)
            fprintf(' Obj   zen[arcsec] azim[deg]  height[m]  lambda[micron] magnitude\n')
            for kObj=1:nObj
                if k>1
                    fprintf(' %2d     %5.2f      %6.2f    [%g %g]     %5.3f          %5.2f\n',...
                        kObj,obj(kObj).zenith*cougarConstants.radian2arcsec,...
                        obj(kObj).azimuth*180/pi,obj(kObj).height,...
                        obj(1,kObj,k).height,...
                        obj(kObj).wavelength*1e6,obj(kObj).magnitude)
                else
                    fprintf(' %2d     %5.2f      %6.2f    %8.2f     %5.3f          %5.2f\n',...
                        kObj,obj(kObj).zenith*cougarConstants.radian2arcsec,...
                        obj(kObj).azimuth*180/pi,obj(kObj).height,...
                        obj(kObj).wavelength*1e6,obj(kObj).magnitude)
                end
            end
            if ~isempty(obj(1).opticalPath)
                fprintf(' Optical path: ')
                cellfun(@(x) fprintf('~>~%s',class(x)), obj(1).opticalPath,'uniformOutput',false);
                fprintf('\n')
            end
            fprintf('----------------------------------------------------\n')
        end
        
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
            obj.directionVector(1) = tan(obj.zenith)*cos(obj.azimuth);
            obj.directionVector(2) = tan(obj.zenith)*sin(obj.azimuth);
            obj.directionVector(3) = 1;
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
        
        function varargout = resetAmplitudeAndPhase(obj)
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
            end
            if nargout>0
                varargout{1} = obj;
            end
        end

        function obj = mtimes(obj,otherObj)
            %% * Source object propagation operator
            %
            % src = src*otherObj propagate src through otherObj
            % multiplying the source amplitude by the otherObj transmitance
            % and adding the otherObj phase to the source phase
            
            nObj = numel(obj);
            for kObj = 1:nObj
                obj(kObj).opticalPath{ length(obj(kObj).opticalPath)+1 } = otherObj;
            end
            relay(otherObj,obj);
        end
        
         function uplus(obj)
            %% UPLUS Source update/stream operator
            %
            % +obj streams the source object through its optical path 
            
            nObj = numel(obj);
            for kObj = 1:nObj
                obj(kObj).mask        = [];
                obj(kObj).p_amplitude = 1;
                obj(kObj).p_phase     = 0;
            end
            cellfun(@(x)relay(x,obj),obj(1).opticalPath,'uniformOutput',false)
        end
              
        function obj = times(obj,otherObj)
            %% .* Source object reset and propagation operator
            %
            % src = src.*otherObj propagate src through otherObj setting the
            % source amplitude to the otherObj transmitance and the source
            % phase to the otherObj phase
            
             mtimes(reset(obj),otherObj);
        end
        
        function out = mldivide(obj,phaseMap)
            %% \ Least square fit to the source wavefront
            %
            % out = obj\phaseMap projects the 2D phase map onto
            % the source wavefront
            
            if isa(phaseMap,'zernike')
                buf = utilities.toggleFrame(obj.phase,2);
                out = buf(obj.mask,:)\phaseMap.p(obj.mask,:);
            else
                out = utilities.toggleFrame(obj.phase,2)\...
                    utilities.toggleFrame(phaseMap,2);
            end
        end

        function out = fresnelPropagation(obj,tel)
            %% FRESNELPROPAGATION Source propagation to the light collector
            %
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
        
        function varargout = polar(obj,lineSpec)
            %% POLAR Display the sources location
            %
            % polar(srcs) plots the sources location around the zenith in
            % arcsec
            %
            % polar(srcs,linespecs) plots the sources location around the
            % zenith in arcsec with given line specification
            %
            % h = polar(...) plots and returns the graphic handle
            
            if nargin<2
                lineSpec = 'o';
            end
            if any(isempty([obj.magnitude]))
                h = polar([obj.azimuth],[obj.zenith]*constants.radian2arcsec,lineSpec);
            else
                h = polar([obj.azimuth],[obj.zenith]*constants.radian2arcsec,'.');
                delete(h)
                nObj = numel(obj);
%                 a = max([obj.magnitude]);
                hold on
                for kObj = 1:nObj
                    h(kObj) = polar(obj(kObj).azimuth,obj(kObj).zenith*constants.radian2arcsec,lineSpec);
                  hcmenu = uicontextmenu;
                  uimenu(hcmenu, 'Label', sprintf('mag=%5.2f',obj(kObj).magnitude))
%                   set(h(kObj),...
%                         'MarkerSize',6*a/obj(kObj).magnitude,...
%                         'MarkerFaceColor',get(h(kObj),'Color'),...
%                         'uicontextmenu',hcmenu);
                end
                hold off
            end
            if nargout>0
                varargout{1} = h;
            end
        end
        
    end
    
    methods (Static)
    
        function [rhoSrcLayer,thetaSrcLayer] = separationMatrix(src1,src2)
            zen1 = [src1.zenith];
            az1  = [src1.azimuth];
            if nargin>1
                zen2 = [src2.zenith];
                az2  = [src2.zenith];
            else
                zen2 = zen1;
                az2 = az1;
            end
            [xSrc1,xSrc2] = meshgrid(tan(zen1).*cos(az1),tan(zen2).*cos(az2));
            [ySrc1,ySrc2] = meshgrid(tan(zen1).*sin(az1),tan(zen2).*sin(az2));
            xSrc = xSrc1 - xSrc2;
            ySrc = ySrc1 - ySrc2;
            rhoSrcLayer = hypot(xSrc,ySrc)*cougarConstants.radian2arcsec;
            thetaSrcLayer = atan2(ySrc,xSrc)*180/pi;
        end
    end
    
end