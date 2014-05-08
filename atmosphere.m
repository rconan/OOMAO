classdef atmosphere < hgsetget
    % Create an atmosphere object
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
    %
    % A full volume of turbulence above a telescope is made combining both
    % a telescope and an atmosphere class
    %
    % The geometric propagation of a wavefront throught the atmosphere is
    % obtained with a source object
    %
    % See also telescope and source
    
    properties
        % Fried parameter
        r0;
        % turbulence outer scale
        L0;
        % number of turbulence layers
        nLayer;
        % turbulence layer object array
        layer;
        % atmosphere tag
        tag = 'ATMOSPHERE';
        % coherence function decay for theta0 and tau0 computation
        % . Roddier (default): coherence function 1/e decay
        % . Fried: structure function equal 1rd^2 (1/\sqrt(e) decay)
        coherenceFunctionDecay = exp(-1);
        % random number generator stream
        rngStream;
        % state at initialization of the atmosphere random number generator stream
        initState;
        wavelengthScale;
        gpu = false;
        % parameters for polar-logarithmic phase screen generation method
        f;
        N_k;
        N_a;
        freq_mag;
        delta_freq_mag;
        freq_ang;
        cos_freq_ang;
        sin_freq_ang;
        sqrt_spectrum_kernel;
        zeta1; 
        eta1;
        zeta2;
        eta2;
    end
    
    properties (Dependent, SetObservable=true)
        % wavelength of the r0
        wavelength;
        % frequency resolution for polar-logarithmic phase screen
        % generation method
        delta;
    end
    
    properties (Dependent, SetAccess=private)
        % seeing
        seeingInArcsec;
        % theta0
        theta0InArcsec;
        % tau0
        tau0InMs;
        % mean height
        meanHeight;
        % mean velocity
        meanWind;
        % Greenwood frequency
        greenwoodFrequency;
    end
    
    properties (Access=private)
        p_wavelength;
        p_log;
        p_delta;
    end
    
    methods
        
        %% Constructor
        function obj = atmosphere(lambda,r0,varargin)
            p = inputParser;
            p.addRequired('wavelength', @(x) isnumeric(x) || isa(x,'photometry') );
            p.addRequired('r0', @isnumeric);
            p.addOptional('L0', Inf, @isnumeric);
            p.addParamValue('altitude', 0, @isnumeric);
            p.addParamValue('fractionnalR0', 1, @isnumeric);
            p.addParamValue('windSpeed', [], @isnumeric);
            p.addParamValue('windDirection', [], @isnumeric);
            p.addParamValue('logging', true, @islogical);
            p.addParamValue('randStream', [], @(x) isempty(x) || isa(x,'RandStream') );
            p.parse(lambda,r0, varargin{:});
            if isa(p.Results.wavelength,'photometry')
                obj.p_wavelength = p.Results.wavelength.wavelength;
            else
                obj.p_wavelength = p.Results.wavelength;
            end
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
            if p.Results.logging
                obj.p_log = logBook.checkIn(obj);
                display(obj);
            end
            if isempty(p.Results.randStream)
                obj.rngStream = RandStream('mt19937ar');
            else
                obj.rngStream = p.Results.randStream;
            end
            obj.initState = obj.rngStream.State;
        end
        
        % Destructor
        function delete(obj)
            if ~isempty(obj.p_log)
                checkOut(obj.p_log,obj)
            end
        end
        
        function newObj = slab(obj,layerIndex)
            %% SLAB Create a single turbulence layer atmosphere object
            %
            % singledAtm = slab(atm,k) creates an atmosphere object from
            % the old atm object and the k-th turbulent layer
            newObj = atmosphere(...
                obj.wavelength,...
                obj.r0,...
                obj.L0,...
                'randStream',obj.rngStream,...
                'logging',false);
            newObj.layer = obj.layer(layerIndex);
        end
        
        function display(obj)
            %% DISPLAY Display object information
            %
            % disp(obj) prints information about the atmosphere object
            
            fprintf('___ %s ___\n',obj.tag)
            if isinf(obj.L0)
                fprintf(' Kolmogorov-Tatarski atmospheric turbulence:\n')
                fprintf('  . wavelength  = %5.2fmicron,\n  . r0          = %5.2fcm,\n  . seeing      = %5.2farcsec,\n  ',...
                    obj.wavelength*1e6,obj.r0*1e2,obj.seeingInArcsec)
            else
                fprintf(' Von Karman atmospheric turbulence\n')
                fprintf('  . wavelength  = %5.2fmicron,\n  . r0          = %5.2fcm,\n  . L0          = %5.2fm,\n  . seeing      = %5.2farcsec,\n  ',...
                    obj.wavelength*1e6,obj.r0*1e2,obj.L0,obj.seeingInArcsec)
            end
            if ~isinf(obj.theta0InArcsec)
                fprintf('. theta0(%2.0f%%) = %5.2farcsec,\n  ',...
                    100*obj.coherenceFunctionDecay,obj.theta0InArcsec)
            end
            if ~isempty(obj.tau0InMs)
                fprintf('. tau0(%2.0f%%)   = %5.2fmillisec',...
                    100*obj.coherenceFunctionDecay,obj.tau0InMs)
            else
%                 fprintf('\b\b')
            end
            fprintf('\n')
            fprintf('----------------------------------------------------\n')
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
            fprintf('----------------------------------------------------\n')
            
        end
        
        %% Set/Get wavelength property
        function val = get.wavelength(obj)
            val = obj.p_wavelength;
        end
        function set.wavelength(obj,val)
            if isa(val,'photometry')
                val = val.wavelength;
            end
            obj.wavelengthScale = obj.wavelength/val;
            obj.r0 = obj.r0.*(val/obj.wavelength)^1.2;
%             if ~isempty(obj.layer(1).phase)
%                 add(obj.p_log,obj,'Scaling wavefront!')
%                 for kLayer=1:obj.nLayer                    
%                     obj.layer(kLayer).phase = ...
%                         obj.layer(kLayer).phase*obj.wavelengthScale;
%                 end
%             end
            obj.p_wavelength = val;
        end
        
        %% Set/Get delta property
        function val = get.delta(obj)
            val = obj.p_delta;
        end
        function set.delta(obj,val)
            l0    = 1e-3;
            L     = 1e2;
            f     = max(L,3*obj.L0);
            kmin  = 2*pi/f;
            f     = f/l0;
            delta = val; % TODO: like f, should depends on the input parameters but 5pct seems optimal from a statistical stand point
            N_k = log(f)/log( (2+delta)/(2-delta) );
            N_a = f.^(1/N_k ) ;
            N_a = 0.25*pi*( (N_a + 1)/(N_a - 1) );
            obj.N_k = ceil( N_k );
            obj.N_a = ceil( N_a );
            obj.f = ( (4*obj.N_a + pi)/(4*obj.N_a - pi) ).^obj.N_k;
            obj.p_delta = 0.5*pi/obj.N_a;
            fprintf('\n@(CEO)>profile:\n');
            fprintf(' . N_k = %f\n',obj.N_k);
            fprintf(' . N_a = %f\n',obj.N_a);
            fprintf(' . frequency range      = %8.2f\n',obj.f);
            fprintf(' . frequency resolution = %6.4f\n',obj.p_delta);
            fprintf(' . minimum frequency    = %6.4f\n',kmin);
            
            ii = (1:obj.N_k)';
            obj.freq_mag = 0.5*kmin*obj.f.^(ii/obj.N_k).*...
                (1 + obj.f.^(-1/obj.N_k) );
            obj.delta_freq_mag = kmin*obj.f.^(ii/obj.N_k).*...
                (1 - obj.f.^(-1/obj.N_k) );
            
            jj = 1:obj.N_a;
            obj.freq_ang = (jj-0.5)*obj.delta;
            obj.cos_freq_ang = cos( obj.freq_ang );
            obj.sin_freq_ang = sin( obj.freq_ang );
            
            k0 = 2*pi/obj.L0;
            obj.sqrt_spectrum_kernel = ...
                ( obj.freq_mag.^2 + k0.^2 ).^(-11/12).*...
                sqrt( obj.freq_mag.*obj.delta_freq_mag.*...
                obj.delta);
            
            obj.zeta1 = rand(obj.N_k,obj.N_a,obj.nLayer);
            obj.zeta1 = sqrt( -log( obj.zeta1 ) );
            obj.eta1  = rand(obj.N_k,obj.N_a,obj.nLayer);
            obj.eta1  = 2.*pi*obj.eta1;
            obj.zeta2 = rand(obj.N_k,obj.N_a,obj.nLayer);
            obj.zeta2 = sqrt( -log( obj.zeta2 ) );
            obj.eta2  = rand(obj.N_k,obj.N_a,obj.nLayer);
            obj.eta2  = 2.*pi*obj.eta2;         
        end
        
        %% Get seeingInArcsec property
        function out = get.seeingInArcsec(obj)
            out = cougarConstants.radian2arcsec.*...
                0.98.*obj.wavelength./obj.r0;
        end
        
        %% Set coherenceFunctionDecay property
        function set.coherenceFunctionDecay(obj,val)
            switch lower(val)
                case 'roddier'
                    obj.coherenceFunctionDecay = exp(-1);
                case 'fried'
                    obj.coherenceFunctionDecay = exp(-0.5);
                otherwise
                    if isnumeric(val)
                        obj.coherenceFunctionDecay = val;
                    else
                        error('oomao:atmosphere:coherenceFunctionDecay',...
                            'The valid definitions for theta0 and tau0 are either Roddier, Fried or a numeric value between 0 and 1 excluded.')
                    end
            end
        end
        
        %% Get theta0InArcsec property
        function out = get.theta0InArcsec(obj)
            z = [obj.layer.altitude];
            if (numel(z)==1 && z==0) || all([obj.layer.altitude]==0)
                out  = Inf;
            else
                if isinf(obj.L0)
                    cst = -log(obj.coherenceFunctionDecay)*(24*gamma(6/5)/5)^(-5/6).*obj.r0.^(5./3);
                    out = ( cst./sum( [obj.layer.fractionnalR0].*z.^(5./3) ) ).^(3/5);
                else
                    fun = @(x) phaseStats.angularStructureFunction(abs(x),obj) + 2.*log(obj.coherenceFunctionDecay);
                    out = abs(fzero(fun,0));
                end
            end
            out = out*cougarConstants.radian2arcsec;
        end
        
        %% Get tau0InMs property
        function out = get.tau0InMs(obj)
            v = [obj.layer.windSpeed];
            if isempty(v)
                out = [];
            elseif numel(v)==1 && v==0
                out  = Inf;
            else
                if isinf(obj.L0)
                    cst = -log(obj.coherenceFunctionDecay)*(24*gamma(6/5)/5)^(-5/6).*obj.r0.^(5./3);
                    out= (cst./sum( [obj.layer.fractionnalR0].*v.^(5./3) ) ).^(3/5);
                else
                    fun = @(x) phaseStats.temporalStructureFunction(abs(x),obj) + 2.*log(obj.coherenceFunctionDecay);
                    out = abs(fzero(fun,0));
                end
            end
            out = out*1e3;
        end
        
        %% Get mean height
        function out = get.meanHeight(obj)
            z = [obj.layer.altitude];
            out = sum( [obj.layer.fractionnalR0].*z.^(5./3) )^(3/5);
        end
        
        %% Get mean wind velocity
        function out = get.meanWind(obj)
            v = [obj.layer.windSpeed];
            out = sum( [obj.layer.fractionnalR0].*v.^(5./3) )^(3/5);
        end
        
        %% Get Greenwood frequency
        function out = get.greenwoodFrequency(obj)
            v = [obj.layer.windSpeed];
            out = 0.4292*...
                sum( [obj.layer.fractionnalR0].*v.^(5./3) )^(3/5)/...
                obj.r0;
        end
        
        function out = polarLogPhaseScreen_(obj,x,y,tau,src)
            %% POLARLOGPHASESCREEN Polar-logarithmic phase screen
            
            if ~any(size(x)==size(y))
                error('OOMAO:atmosphere:polarLogPhaseScreen',...
                    'x and y must have the same size!');
            end
            n = numel(x); 
            hl = zeros(1,1,obj.nLayer);
            vx = zeros(1,1,obj.nLayer);
            vy = zeros(1,1,obj.nLayer);
            xi0 = zeros(1,1,obj.nLayer);
            hl(1,1,:) = [obj.layer.altitude];
            [vx(1,1,:),vy(1,1,:)] = pol2cart( ...
                [obj.layer.windDirection],...
                [obj.layer.windSpeed]);
            xi0(1,1,:) = sqrt( [obj.layer.fractionnalR0] );
            gl = 1 - hl./src.height;
            directionVector = src.directionVector;
            out = zeros( size(x) );
            parfor k=1:n
                xl = gl.* ( x(k) - vx*tau ) + ...
                    hl*directionVector(1);
                yl = gl.* ( y(k) - vy*tau ) + ...
                    hl*directionVector(2);
                red = obj.zeta1.*cos( bsxfun( @plus, obj.eta1, ...
                    bsxfun( @times, obj.freq_mag, ...
                    bsxfun( @times, xl , obj.cos_freq_ang) + ...
                    bsxfun( @times, yl , obj.sin_freq_ang) ) ) ) + ...
                    obj.zeta2.*cos( bsxfun( @minus, obj.eta2, ...
                    bsxfun( @times, obj.freq_mag, ...
                    bsxfun( @times, xl , obj.sin_freq_ang) - ...
                    bsxfun( @times, yl , obj.cos_freq_ang) ) ) );
                red = bsxfun( @times, xi0, red);
                red = sum(red,3);
                red = bsxfun( @times, ...
                    obj.sqrt_spectrum_kernel, red);
                out(k) = sum(red(:));
            end
            out = 1.4.*obj.r0.^(-5/6).*out;
        end
        
        function out = fourierPhaseScreen(atm,D,nPixel,nMap)
            %% FOURIERPHASESCREEN Phase screen computation
            %
            % map = fourierPhaseScreen(atm,D,nPixel) Computes a square
            % phase screen of D meter and sampled with nPixel using the
            % Fourier method
            %
            % map = fourierPhaseScreen(atm) Computes a square phase screen
            % of atm.layer.D meter and sampled with atm.layer.nPixel using
            % the Fourier method; atm must contain only one turbulent
            % layer!
            %
            % See also atmosphere
            
            %             warning('oomao:atmosphere:fourierPhaseScreen',...
            %                 'The fourierPhaseScreen seems to have a bug, to use with care!')
            
            if nargin<4
                nMap = 1;
            end
            if nargin<2
                D = atm.layer.D;
                nPixel = atm.layer.nPixel;
            end
            
            N = 4*nPixel;
            L = (N-1)*D/(nPixel-1);
            [fx,fy]  = freqspace(N,'meshgrid');
            [~,fr]  = cart2pol(fx,fy);
            fr  = fftshift(fr.*(N-1)/L./2);
            clear fx fy fo
            psdRoot = sqrt(phaseStats.spectrum(fr,atm)); % Phase FT magnitude
            %             figure
            %             imagesc(map)
            %             axis square
            %             colorbar
            clear fr
            fourierSampling = 1./L;
            %                         % -- Checking the variances --
            %                         theoreticalVar = phaseStats.variance(atm);
            %                         disp(['Info.: Theoretical variance: ',num2str(theoreticalVar,'%3.3f'),'rd^2'])
            %                         %     numericalVar    = sum(abs(phMagSpectrum(:)).^2).*fourierSampling.^2;
            %                         numericalVar    = trapz(trapz(map.^2)).*fourierSampling.^2;
            %                         disp(['Info.: Numerical variance  :',num2str(numericalVar,'%3.3f'),'rd^2'])
            %                         % -------------------------------
            u = 1:nPixel;
            out = zeros(nPixel,nPixel,nMap);
            for kMap=1:nMap
                map = psdRoot.*fft2(randn(atm.rngStream,N))./N; % White noise filtering
                map = real(ifft2(map).*fourierSampling).*N.^2;
                out(:,:,kMap) = map(u,u);
            end
        end
        function out = fourierPhaseScreenStraight(atm,D,nPixel)
            %% FOURIERPHASESCREEN Phase screen computation
            %
            % map = fourierPhaseScreen(atm,D,nPixel) Computes a square
            % phase screen of D meter and sampled with nPixel using the
            % Fourier method
            %
            % map = fourierPhaseScreen(atm) Computes a square phase screen
            % of atm.layer.D meter and sampled with atm.layer.nPixel using
            % the Fourier method; atm must contain only one turbulent
            % layer!
            %
            % See also atmosphere
            
            %             warning('oomao:atmosphere:fourierPhaseScreen',...
            %                 'The fourierPhaseScreen seems to have a bug, to use with care!')
            
            if nargin<2
                D = atm.layer.D;
                nPixel = atm.layer.nPixel;
            end
            
            N = nPixel;
            del_f = 1/D;   % frequency grid spacing [1/m]
            fx = (-N/2 : N/2-1) * del_f;
            % frequency grid [1/m]
            [fx fy] = meshgrid(fx);
            [~, f] = cart2pol(fx, fy);  % polar grid
            % fm = 5.92/l0/(2*pi); % inner scale frequency [1/m]
            % f0 = 1/L0;           % outer scale frequency [1/m]
            % von Karman atmospheric phase PSD
            % PSD_phi = 0.023*r0^(-5/3) * exp(-(f/fm).^2) ...
            %     ./ (f.^2 + f0^2).^(11/6);
            PSD_phi = phaseStats.spectrum(f,atm);
            PSD_phi(N/2+1,N/2+1) = 0;
            % random draws of Fourier coefficients
            cn = (randn(atm.rngStream,N) + 1i*randn(atm.rngStream,N)) .* sqrt(PSD_phi)*del_f;
            % synthesize the phase screen
            out = real(ifftshift(ifft2(ifftshift(cn))))*N^2;
            
        end
        
        function out = fourierSubHarmonicPhaseScreen(atm,D,nPixel,nMap)
            %% FOURIERSUBHARMONICPHASESCREEN Phase screen computation
            %
            % map = fourierPhaseScreen(atm,D,nPixel) Computes a square
            % phase screen of D meter and sampled with nPixel using the
            % Fourier method
            %
            % map = fourierPhaseScreen(atm) Computes a square phase screen
            % of atm.layer.D meter and sampled with atm.layer.nPixel using
            % the Fourier method; atm must contain only one turbulent
            % layer!
            %
            % See also atmosphere
            
            %             warning('oomao:atmosphere:fourierPhaseScreen',...
            %                 'The fourierPhaseScreen seems to have a bug, to use with care!')
            
            if nargin>3
                out = zeros(nPixel,nPixel,nMap);
                h = waitbar(0,'Phase screens computing ...');
                for kMap=1:nMap
                    out(:,:,kMap) = fourierSubHarmonicPhaseScreen(atm,D,nPixel);
                    waitbar(kMap/nMap,h)
                end
                close(h)
            else
            if nargin<2
                D = atm.layer.D;
                nPixel = atm.layer.nPixel;
            end
            
            % function [phz_lo phz_hi] ...
            %     = ft_sh_phase_screen(r0, N, delta, L0, l0)
            
            %     D = N*delta;
            N = nPixel;
            delta = D/N;
            % high-frequency screen from FFT method
            %     phz_hi = ft_phase_screen(atm.r0, N, delta, atm.L0, 0);
            phz_hi = fourierPhaseScreenStraight(atm,D,nPixel);
            % spatial grid [m]
            [x y] = meshgrid((-N/2 : N/2-1) * delta);
            % initialize low-freq screen
            phz_lo = zeros(size(phz_hi));
            % loop over frequency grids with spacing 1/(3^p*L)
            for p = 1:3
                % setup the PSD
                del_f = 1 / (3^p*D); % frequency grid spacing [1/m]
                fx = (-1 : 1) * del_f;
                % frequency grid [1/m]
                [fx fy] = meshgrid(fx);
                [~, f] = cart2pol(fx, fy);  % polar grid
%                 fm = 5.92/l0/(2*pi); % inner scale frequency [1/m]
%                 f0 = 1/L0;           % outer scale frequency [1/m]
                % modified von Karman atmospheric phase PSD
                %         PSD_phi = 0.023*r0^(-5/3) * exp(-(f/fm).^2) ...
                %             ./ (f.^2 + f0^2).^(11/6);
                PSD_phi = phaseStats.spectrum(f,atm);
                PSD_phi(2,2) = 0;
                % random draws of Fourier coefficients
                cn = (randn(atm.rngStream,3) + 1i*randn(atm.rngStream,3)) ...
                    .* sqrt(PSD_phi)*del_f;
                SH = zeros(N);
                % loop over frequencies on this grid
                for ii = 1:9
                    SH = SH + cn(ii) ...
                        * exp(i*2*pi*(fx(ii)*x+fy(ii)*y));
                end
                phz_lo = phz_lo + SH;   % accumulate subharmonics
            end
            phz_lo = real(phz_lo) - mean(real(phz_lo(:)));
            out = phz_hi + phz_lo;
            end
        end
        
        function map = choleskyPhaseScreen(atm,D,nPixel,nMap)
            %% CHOLESKYPHASESCREEN Phase screen computation 
            %
            % map = choleskyPhaseScreen(obj,D,nPixel) Computes a square
            % phase screen of D meter and sampled with nPixel using the
            % Cholesky factorisation of the atmospheric phase covariance
            % matrix
            %
            % map = choleskyPhaseScreen(obj,D,nPixel,covarianceMatrix)
            % Computes a square phase screen of D meter and sampled with
            % nPixel using the Cholesky factorisation of the atmospheric
            % phase covariance matrix covarianceMatrix
            %
            % See also chol and atmosphere
            
%             if nargin<2
%                 D = atm.layer.D;
%                 nPixel = atm.layer.nPixel;
%             end
            if nargin<4
                nMap = 1;
            end
            if nargin==1
                for kLayer=1:atm.nLayer
                    
                    nPixel = atm.layer(kLayer).nPixel;
                    if isempty(atm.layer(kLayer).choleskyFact)
                        fprintf('Computing the Cholesky factor matrix!\n')
                        D = atm.layer(kLayer).D;
%                         [x,y] = meshgrid((0:nPixel-1)*D/nPixel);
                        [x,y] = meshgrid(linspace(-1,1,nPixel)*D/2);
                        atm.layer(kLayer).choleskyFact = chol( ...
                            phaseStats.covarianceToeplitzMatrix(slab(atm,kLayer),complex(x,y)) ,'lower');
                    end
                    map = atm.layer(kLayer).choleskyFact*randn(atm.rngStream,nPixel^2,nMap);
                    map = reshape(map,nPixel,nPixel,nMap);
                    atm.layer(kLayer).phase = map;
                    
                end
            else
                fprintf('Computing the Cholesky factor matrix!\n')
%                 [x,y] = meshgrid((0:nPixel-1)*D/nPixel);
                [x,y] = meshgrid(linspace(-1,1,nPixel)*D/2);
                choleskyFact = chol( ...
                    phaseStats.covarianceToeplitzMatrix(atm,complex(x,y)) ,'lower');
                map = choleskyFact*randn(atm.rngStream,nPixel^2,nMap);
                map = reshape(map,nPixel,nPixel,nMap);
            end
        end
        
        function host2gpu(obj)
            host2gpu(obj.layer)
            obj.gpu = true;
        end
        
        function gpu2host(obj)
            gpu2host(obj.layer)
            obj.gpu = false;
        end
        
    end
    
    methods (Static)
        
        function out = r0VsWavelength(r0Wavelength,r0,wavelength)
            %% R0VSWAVELENGTH r0 versus wavelength
            %
            % out = r0VsWavelength(r0Wavelength,r0,wavelength) computes the
            % value of r0(r0Wavelength) at the new wavelength
        
            if ischar(r0Wavelength)
                r0Wavelength = eval(['photometry.',r0Wavelength]);
            end
            if ischar(wavelength)
                wavelength = eval(['photometry.',wavelength]);
            end
            out = r0.*(wavelength./r0Wavelength).^1.2;
            
        end
        
    end
    
end

function phz = ft_phase_screen(r0, N, delta, L0, l0)
% function phz ...
%     = ft_phase_screen(r0, N, delta, L0, l0)

% setup the PSD
del_f = 1/(N*delta);   % frequency grid spacing [1/m]
fx = (-N/2 : N/2-1) * del_f;
% frequency grid [1/m]
[fx fy] = meshgrid(fx);
[~, f] = cart2pol(fx, fy);  % polar grid
fm = 5.92/l0/(2*pi); % inner scale frequency [1/m]
f0 = 1/L0;           % outer scale frequency [1/m]
% modified von Karman atmospheric phase PSD
PSD_phi = 0.023*r0^(-5/3) * exp(-(f/fm).^2) ...
    ./ (f.^2 + f0^2).^(11/6);
PSD_phi(N/2+1,N/2+1) = 0;
% random draws of Fourier coefficients
cn = (randn(N) + 1i*randn(N)) .* sqrt(PSD_phi)*del_f;
% synthesize the phase screen
phz = real(ifftshift(ifft2(ifftshift(cn))))*N^2;
end
