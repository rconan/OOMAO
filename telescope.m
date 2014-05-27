classdef telescope < telescopeAbstract
    % Create a telescope object
    %
    % sys = telescope(D) creates a telescope object from the telescope diameter D
    %
    % sys = telescope(D,...) creates a telescope object from the
    % above parameter and from optionnal parameter-value pair arguments. The
    % optionnal parameters are those of the telescopeAbstract class.
    %
    % Example:
    % tel = telescope(8); An 8m diameter telescope
    %
    % tel = telescope(8,'obstructionRatio',0.14,'fieldOfViewInArcmin',2,'resolution',64); 
    % An 8m diameter telescope with an 14% central obstruction, a 2 arcmin
    % fov and a pupil sampled with 64x64 pixels
    % Displaying the pupil:
    % imagesc(tel.pupil)
    % The pupil logical mask is given by
    % tel.pupilLogical
    %
    % A telescope object can be combined with an atmosphere object to
    % define a volume of turbulence above the telescope within its
    % field-of--view:
    % atm = atmosphere(photometry.V,0.15,30,...
    %     'altitude',4e3,...
    %     'fractionnalR0',1,...
    %     'windSpeed',15,...
    %     'windDirection',0);
    % tel = telescope(8,'fieldOfViewInArcmin',2,'resolution',64,'samplingTime',1/500); 
    % The telescope-atmosphere object is built by adding the atmosphere to
    % the telescope:
    % tel = tel + atm; 
    % figure, imagesc(tel)
    % The frozen-flow or Taylor motion of the phase screen is created by
    % updating the telescope object:
    % +tel;
    % The phase screen(s) are moved of an amount depending on the
    % samplingTime and the wind vector parameters.
    % The geometric propagation of a star through the atmosphere to the
    % telescope pupil is done as followed
    % ngs = source;
    % ngs = ngs.*tel;
    % figure, imagesc(ngs.phase), axis square xy, colorbar
    % At anytime the atmosphere can be removed from the telescope:
    % tel = tel - atm;
    %
    % See also telescopeAbstract, atmosphere and source
    
    properties
        % telescope tag
        tag = 'TELESCOPE';
    end
    
    
    properties (Dependent)% , SetAccess = private)
        % telescope pupil mask
        pupil;
    end
    
%     properties (Access=private)
%         atm;
%         innerMask;
%         outerMask;
%         A;
%         B;
%         windVx;
%         windVy;
%         count;
%         mapShift;
%         nShift;
%         x;
%         y;
%         imageHandle;
%         layerSampling;
%         sampler;
%         log;
%         p_pupil;
%     end
    
    methods
        
        %% Constructor
        function obj = telescope(D,varargin)
%             p = inputParser;
%             p.addRequired('D', @isnumeric);
%             p.addParamValue('obstructionRatio', 0, @isnumeric);
%             p.addParamValue('fieldOfViewInArcsec', [], @isnumeric);
%             p.addParamValue('fieldOfViewInArcmin', [], @isnumeric);
%             p.addParamValue('resolution', [], @isnumeric);
%             p.addParamValue('samplingTime', [], @isnumeric);
%             p.addParamValue('opticalAberration', [], @(x) isa(x,'atmosphere'));
%             p.parse(D,varargin{:});
            obj = obj@telescopeAbstract(D,varargin{:});
            display(obj)
        end
        
        %% Destructor
        function delete(obj)
            if isa(obj.opticalAberration,'atmosphere')
                add(obj.log,obj,'Deleting atmosphere layer slabs!')
%                 for kLayer=1:obj.atm.nLayer
%                     obj.atm.layer(kLayer).phase = [];
%                 end
            end
            checkOut(obj.log,obj)
        end
        
        function obj = saveobj(obj)
            obj.phaseListener = [];
        end
        
        function display(obj)
            %% DISPLAY Display object information
            %
            % display(obj) prints information about the atmosphere+telescope object
            
            display(obj.atm)
            fprintf('___ %s ___\n',obj.tag)
            if obj.obstructionRatio==0
                fprintf(' %4.2fm diameter full aperture',obj.D)
            else
                fprintf(' %4.2fm diameter with a %4.2f%% central obstruction',...
                    obj.D,obj.obstructionRatio*100)
            end
            fprintf(' with %5.2fm^2 of light collecting area;\n',obj.area)
            if obj.fieldOfView~=0
                fprintf(' the field-of-view is %4.2farcmin;',...
                    obj.fieldOfView*constants.radian2arcmin)
            end
            if ~isempty(obj.resolution)
                fprintf(' the pupil is sampled with %dX%d pixels',...
                    obj.resolution,obj.resolution)
            end
            if obj.fieldOfView~=0 || ~isempty(obj.resolution)
                fprintf('\n')
            end
            fprintf('----------------------------------------------------\n')
            
        end

        %% Get and Set the pupil
        function pupil = get.pupil(obj)
            pupil = obj.p_pupil;
            if isempty(pupil) && ~isempty(obj.resolution)
                pupil = utilities.piston(obj.resolution,'shape',obj.shape);
                if obj.obstructionRatio>0
                    pupil = pupil - ...
                        utilities.piston(...
                        round(obj.resolution.*obj.obstructionRatio),...
                        obj.resolution);
                end
                obj.p_pupil = pupil;
            end
        end
        function obj = set.pupil(obj,val)
            obj.p_pupil = val;
        end
        
                 
        function out = otf(obj, r)
            %% OTF Telescope optical transfert function
            %
            % out = otf(obj, r) Computes the telescope optical transfert function
            
            if ~all(isreal(r(:))) % Check if r is complex!
                r = abs(r);
            end
%             out = zeros(size(r));
            if obj.obstructionRatio ~= 0
                out = pupAutoCorr(obj.D) + pupAutoCorr(obj.obstructionRatio*obj.D) - ...
                    2.*pupCrossCorr(obj.D./2,obj.obstructionRatio*obj.D./2);
            else
                out = pupAutoCorr(obj.D);
            end
            out = out./(pi*obj.D*obj.D.*(1-obj.obstructionRatio*obj.obstructionRatio)./4);
            
            if isa(obj.opticalAberration,'atmosphere')
                out = out.*phaseStats.otf(r,obj.opticalAberration);
            end
            
            function out1 = pupAutoCorr(D)
                
                index       = r <= D;
                red         = r(index)./D;
                out1        = zeros(size(r));
                out1(index) = D.*D.*(acos(red)-red.*sqrt((1-red.*red)))./2;
                
            end            
            
            function out2 = pupCrossCorr(R1,R2)
                
                out2 = zeros(size(r));
                
                index       = r <= abs(R1-R2);
                out2(index) = pi*min([R1,R2]).^2;
                
                index       = (r > abs(R1-R2)) & (r < (R1+R2));
                rho         = r(index);
                red         = (R1*R1-R2*R2+rho.*rho)./(2.*rho)/(R1);
                out2(index) = out2(index) + R1.*R1.*(acos(red)-red.*sqrt((1-red.*red)));
                red         = (R2*R2-R1*R1+rho.*rho)./(2.*rho)/(R2);
                out2(index) = out2(index) + R2.*R2.*(acos(red)-red.*sqrt((1-red.*red)));
                
            end
            
        end
         
        function out = psf(obj,f,~)
            %% PSF Telescope point spread function
            %
            % out = psf(obj, f) computes the telescope point spread function
            
            if isa(obj.opticalAberration,'atmosphere')
                fun = @(u) 2.*pi.*quadgk(@(v) psfHankelIntegrandNested(v,u),0,obj.D);
                out = zeros(size(f));
                parfor k = 1:numel(f)
                    out(k) = fun(f(k));
                end
%                 out = arrayfun( fun, f);
            else
                out   = ones(size(f)).*pi.*obj.D.^2.*(1-obj.obstructionRatio.^2)./4;
                index = f~=0;
                u = pi.*obj.D.*f(index);
                surface = pi.*obj.D.^2./4;
                out(index) = surface.*2.*besselj(1,u)./u;
                if obj.obstructionRatio>0
                    u = pi.*obj.D.*obj.obstructionRatio.*f(index);
                    surface = surface.*obj.obstructionRatio.^2;
                    out(index) = out(index) - surface.*2.*besselj(1,u)./u;
                end
                out = abs(out).^2./(pi.*obj.D.^2.*(1-obj.obstructionRatio.^2)./4);
                
            end
            function y = psfHankelIntegrandNested(x,freq)
                y = x.*besselj(0,2.*pi.*x.*freq).*otf(obj,x);
            end
        end
        
        function out = image(obj,resolution,pixelScaleInSpFreq)
            %% IMAGE 2D Point Spread Function
            %
            % psf = image(tel,resolution,pixelScaleInSpFreq)
            
            n = resolution;
            u = pixelScaleInSpFreq*linspace(-1,1,n)*n/2;
            [fx,fy] = meshgrid(u);
            out = psf(obj,hypot(fx,fy))/obj.area;%/psf(obj,0);

        end
             
        function out = fullWidthHalfMax(obj)
            %% FULLWIDTHHALFMAX Full Width at Half the Maximum evaluation
            %
            % out = fullWidthHalfMax(a) computes the FWHM of a telescope
            % object. Units are m^{-1}. To convert it in arcsecond,
            % multiply by the wavelength then by radian2arcsec.
            
            if isa(obj.opticalAberration,'atmosphere')
                x0 = [0,2/min(obj.D,obj.opticalAberration.r0)];
            else
                x0 = [0,2/obj.D];
            end
            [out,fval,exitflag] = fzero(@(x) psf(obj,abs(x)./2) - psf(obj,0)./2,x0,optimset('TolX',1e-9));
            if exitflag<0
                warning('cougar:telescope:fullWidthHalfMax',...
                    'No interval was found with a sign change, or a NaN or Inf function value was encountered during search for an interval containing a sign change, or a complex function value was encountered during the search for an interval containing a sign change.')
            end
            out = abs(out);
        end
        
        function varargout = footprintProjection(obj,zernModeMax,src)
            nSource = length(src);
            P = cell(obj.atm.nLayer,nSource);
            obj.log.verbose = false;
            for kSource = 1:nSource
                fprintf(' @(telescope) > Source #%2d - Layer #00',kSource)
                for kLayer = 1:obj.atm.nLayer
                    fprintf('\b\b%2d',kLayer)
                    obj.atm.layer(kLayer).zern = ...
                        zernike(1:zernModeMax,'resolution',obj.atm.layer(kLayer).nPixel);
                    conjD = obj.atm.layer(kLayer).D;
                    delta = obj.atm.layer(kLayer).altitude.*...
                        tan(src(kSource).zenith).*...
                        [cos(src(kSource).azimuth),sin(src(kSource).azimuth)];
                    delta = delta*2/conjD;
                    alpha = conjD./obj.D;
                    P{kLayer,kSource} = smallFootprintExpansion(obj.atm.layer(kLayer).zern,delta,alpha);
                    varargout{1} = P;
                end
                fprintf('\n')
            end
            obj.log.verbose = true;
%             if nargout>1
%                 o = linspace(0,2*pi,101);
%                 varargout{2} = cos(o)./alpha + delta(1);
%                 varargout{3} = sin(o)./alpha + delta(2);
%             end
        end
        
        function out = index(obj,z)
            %% INDEX
            
            out = abs(z)<obj.R & abs(z)>obj.R*obj.obstructionRatio;
        end
        
    end
    
        
    methods (Static)
        function sys = demo(action)
            sys = telescope(2.2e-6,0.8,30,25,...
                'altitude',[0,10,15].*1e3,...
                'fractionnalR0',[0.7,0.2,0.1],...
                'windSpeed',[10,5,15],...
                'windDirection',[0,pi/4,pi/2],...
                'fieldOfViewInArcMin',2,...
                'resolution',60*15,...
                'samplingTime',1/500);
            %             sys = telescope(2.2e-6,0.8,10,25,...
            %                 'altitude',10.*1e3,...
            %                 'fractionnalR0',1,...
            %                 'windSpeed',10,...
            %                 'windDirection',0,...
            %                 'fieldOfViewInArcMin',1,...
            %                 'resolution',60,...
            %                 'samplingTime',1e-3);
            if nargin>0
                update(sys);
                sys.phaseListener.Enabled = true;
                imagesc(sys)
                while true
                    update(sys);
                    drawnow
                end
            end
        end
        
        function sys = demoSingleLayer
            sys = telescope(2.2e-6,0.8,10,25,...
                'altitude',10e3,...
                'fractionnalR0',1,...
                'windSpeed',10,...
                'windDirection',0,...
                'fieldOfViewInArcMin',2,...
                'resolution',256,...
                'samplingTime',1/500);
        end
    end
    
end
