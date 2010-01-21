classdef atmosphere < handle
% ATMOSPHERE Create an atmosphere object
%
% atm = atmosphere(wavelength,r0) creates an atmosphere object from the
% wavelength and the Fried parameter r0
%
% atm = atmosphere(wavelength,r0,L0) creates an atmosphere object from the
% wavelength, the Fried parameter r0 and the outer scale L0
%
% atm =
% atmosphere(wavelength,r0,'altitude',altitude,'fractionnalR0',fractionnalR
% 0) creates an atmosphere object from the wavelength, the Fried parameter
% r0, and from the altitudes and the fractionnalR0 of the turbulence layers
%
% atm =
% atmosphere(wavelength,r0,'altitude',altitude,'fractionnalR0',fractionnalR
% 0,'windSpeed',windSpeed,'windDirection',windDirection) creates an
% atmosphere object from the wavelength, the Fried parameter r0, and from
% the altitudes, the fractionnalR0, the wind speeds and the wind directions
% of the turbulence layers
%
% atm =
% atmosphere(wavelength,r0,L0,'altitude',altitude,'fractionnalR0',fractionn
% alR0) creates an atmosphere object from the wavelength, the Fried
% parameter r0, the outer scale L0 and from the altitudes and the
% fractionnalR0 of the turbulence layers
%
% atm =
% atmosphere(wavelength,r0,L0,'altitude',altitude,'fractionnalR0',fractionn
% alR0,'windSpeed',windSpeed,'windDirection',windDirection) creates an
% atmosphere object from the wavelength, the Fried parameter r0, the outer
% scale L0 and from the altitudes, the fractionnalR0, the wind speeds and
% the wind directions of the turbulence layers
%
% Example:
%     atm = atmosphere(photometry.V,0.15,30,...
%     'altitude',4e3,...
%     'fractionnalR0',1,...
%     'windSpeed',15,...
%     'windDirection',0);
% atm.wavelength = photometry.R;

    properties
        % Fried parameter
        r0;
        % turbulence outer scale
        L0;
        % number of turbulence layers
        nLayer;
        % turbulence layer object array
        layer;
    end
    
    properties (Dependent)
        % wavelength of the r0
        wavelength;
    end
    
    properties (Dependent, SetAccess=private)
        % seeing
        seeingInArcsec;
        % theta0
        theta0InArcsec;
        % tau0
        tau0InMs;
    end
    
    properties (Access=private)
        p_wavelength;
%         log;
    end
    
    methods
        
        %% Constructor
        function obj = atmosphere(lambda,r0,varargin)
            p = inputParser;
            p.addRequired('wavelength', @isnumeric);
            p.addRequired('r0', @isnumeric);
            p.addOptional('L0', Inf, @isnumeric);
            p.addParamValue('altitude', 0, @isnumeric);
            p.addParamValue('fractionnalR0', 1, @isnumeric);
            p.addParamValue('windSpeed', [], @isnumeric);
            p.addParamValue('windDirection', [], @isnumeric);
            p.parse(lambda,r0, varargin{:});
            obj.p_wavelength = p.Results.wavelength;
            obj.r0 = p.Results.r0;
            obj.L0 = p.Results.L0;
            obj.nLayer = length(p.Results.altitude);
            if any(isempty(p.Results.windSpeed))
                obj.layer = turbulenceLayer(...
                    p.Results.altitude,...
                    p.Results.fractionnalR0);
            else
                obj.layer = turbulenceLayer(...
                    p.Results.altitude,...
                    p.Results.fractionnalR0,...
                    p.Results.windSpeed,...
                    p.Results.windDirection);
            end
%             obj.log = logBook.checkIn(obj);
        end
        
%         % Destructor
%         function delete(obj)
%             checkOut(obj.log,obj)
%         end
        
        function newObj = slab(obj,layerIndex)
            %% SLAB Create a single turbulence layer atmosphere object
            %
            % singledAtm = slab(atm,k) creates an atmosphere object from
            % the old atm object and the k-th turbulent layer
            newObj = atmosphere(...
                obj.wavelength,...
                obj.r0,...
                obj.L0);
            newObj.layer = obj.layer(layerIndex);
        end
        
        function display(obj)
            %% DISP Display object information
            %
            % disp(obj) prints information about the atmosphere object
            
            if isinf(obj.L0)
                fprintf(' Kolmogorov-Tatarski atmospheric turbulence:\n')
                fprintf('  wavelength=%.2fmicron, r0=%.2fcm, seeing=%.2farcsec,\n  ',...
                obj.wavelength*1e6,obj.r0*1e2,obj.seeingInArcsec)
            else
                fprintf(' Von Karman atmospheric turbulence\n')
                fprintf('  wavelength=%.2fmicron, r0=%.2fcm, L0=%.2fm, seeing=%.2farcsec,\n  ',...
                obj.wavelength*1e6,obj.r0*1e2,obj.L0,obj.seeingInArcsec)
            end
            if ~isinf(obj.theta0InArcsec) 
                fprintf('theta0=%.2farcsec, ',obj.theta0InArcsec)
            end
            if ~isempty(obj.tau0InMs)
                 fprintf('tau0=%.2fmillisec',obj.tau0InMs)
            else
                fprintf('\b\b')
            end
            fprintf('\n')
                fprintf('  Layer   Altitude[m]   fr0    wind([m/s] [deg])   D[m]    res[px]\n')
            for kLayer=1:obj.nLayer
                fprintf('  %2d      %8.2f      %4.2f    (%5.2f %6.2f)     %5.2f    %3d\n',...
                    kLayer,...
                    obj.layer(kLayer).altitude,...
                    obj.layer(kLayer).fractionnalR0,...
                    obj.layer(kLayer).windSpeed,...
                    obj.layer(kLayer).windDirection*180/pi,...
                    obj.layer(kLayer).D,...
                    obj.layer(kLayer).nPixel)
            end
       end
        function val = get.wavelength(obj)
            val = obj.p_wavelength;
        end
        
        function set.wavelength(obj,val)
            obj.r0 = obj.r0.*(val/obj.wavelength)^1.2;
            obj.p_wavelength = val;
        end
        
        function out = get.seeingInArcsec(obj)
            out = cougarConstants.radian2arcsec.*...
                0.98.*obj.wavelength./obj.r0;
        end
        
        function out = get.theta0InArcsec(obj)
            z = [obj.layer.altitude];
            if numel(z)==1 && z==0
               out  = Inf;
            else
                if isinf(obj.L0)
                    cst = (24*gamma(6/5)/5)^(-5/6).*obj.r0.^(5./3);
                    out= ( cst./sum( [obj.layer.fractionnalR0].*z.^(5./3) ) ).^(3/5);
                else
                    fun = @(x) phaseStats.angularStructureFunction(abs(x),obj) + 2.*log(exp(-1));
                    out = abs(fzero(fun,0));
                end
            end
            out = out*cougarConstants.radian2arcsec;
        end
        
        function out = get.tau0InMs(obj)
            v = [obj.layer.windSpeed];
            if isempty(v)
                out = [];
            elseif numel(v)==1 && v==0
               out  = Inf;
            else
                if isinf(obj.L0)
                    cst = (24*gamma(6/5)/5)^(-5/6).*obj.r0.^(5./3);
                    out= ( cst./sum( [obj.layer.fractionnalR0].*v.^(5./3) ) ).^(3/5);
                else
                    fun = @(x) phaseStats.temporalStructureFunction(abs(x),obj) + 2.*log(exp(-1));
                    out = abs(fzero(fun,0));
                end
            end
            out = out*1e3;
        end
        
        function map = fourierPhaseScreen(atm,D,nPixel)
            N = 4*nPixel;
            L = (N-1)*D/(nPixel-1);
            [fx,fy]  = freqspace(N,'meshgrid');
            [fo,fr]  = cart2pol(fx,fy);
            fr  = fftshift(fr.*(N-1)/L./2);
            clear fx fy fo
            map = sqrt(phaseStats.spectrum(fr,atm)); % Phase FT magnitude
            clear fr
            fourierSampling = 1./L;
            %             % -- Checking the variances --
            %             theoreticalVar = variance(phaseStats1);
            %             disp(['Info.: Theoretical variance: ',num2str(theoreticalVar,'%3.3f'),'rd^2'])
            %             %     numericalVar    = sum(abs(phMagSpectrum(:)).^2).*fourierSampling.^2;
            %             numericalVar    = trapz(trapz(amp.^2)).*fourierSampling.^2;
            %             disp(['Info.: Numerical variance  :',num2str(numericalVar,'%3.3f'),'rd^2'])
            %             % -------------------------------
            map = map.*fft2(randn(N))./N; % White noise filtering
            map = real(ifft2(map).*fourierSampling).*N.^2;
            u = 1:nPixel;
            map = map(u,u);
        end
        
    end
    
end