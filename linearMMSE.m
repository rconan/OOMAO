classdef linearMMSE < handle
    %% LINEARMMSE Create a linearMMSE object
    %
    % mmse = linearMMSE(sampling,diameter,atmModel,guideStar)
    % mmse = linearMMSE(sampling,diameter,atmModel,guideStar,mmseStar)
    % mmse = linearMMSE(...,'pupil',pupilMask,'unit',unit)
    % Results are given at the wavelength of the mmseStar
    %
    % Example:
    % sampling = 51;
    % diameter = 25.4;
    % atmModel = gmtAtmosphere(1,60);
    % gs = source('asterism',{[6,arcsec(35),0]},'wavelength',photometry.Na);
    % ss = source('wavelength',photometry.H);
    % mmse = linearMMSE(sampling,diameter,atmModel,gs,ss,'pupil',utilities.piston(51,'type','logical'));
    
    properties
        tag = 'linearMMSE';
        % pupil sampling
        sampling
        % pupil diameter
        diameter
        % pupil mask
        pupil
        % atmosphere mode
        atmModel
        % guide stars
        guideStarListener
        % guide star covariance matrix
        Cxx
        % mmse star covariance matrix
        Coo
        % mmse star / guide star covariance matrix
        Cox
        nGuideStar
        nmmseStar
        % error unit
        unit;
        % Bayesian minimum mean square error        
        Bmse;
        % zonal of modal model
        model;
        % Zernike mode for zonal model
        zernikeMode;
        % zernike / guide star covariance matrix
        Cax
        % zernike / mmse star covariance matrix
        Cao
        % zernike polynomials matrix
        Z
        % zernike covariance matrix
        Caa
        % mmse reconstructor
        mmseBuilder
        % Zernike noise variance
        zernNoiseCovariance
    end
    
    properties (SetObservable=true)
        guideStar
    end
    
    properties (Dependent)
        % mmse estimates directions
        mmseStar
    end
    
    properties (Dependent,SetAccess=private)
        % error variance
        var
        % error rms
        rms
        % error rms in mas
        rmsMas
        % error variance map
        varMap
        % error rms map
        rmsMap        
    end
    
    properties (Dependent)
        %  noise covariance
        noiseCovariance
    end
    
    properties (Access=private)
        log
        p_mmseStar
        zernP
        p_noiseCovariance = 0;
        tel;
    end
 
    methods
        
        %% Constructor
        function obj = linearMMSE(sampling,diameter,atmModel,guideStar,varargin)
            
            inputs = inputParser;
            inputs.addRequired('sampling',@isnumeric);
            inputs.addRequired('diameter',@isnumeric);
            inputs.addRequired('atmModel',@(x) isa(x,'atmosphere'));
            inputs.addRequired('guideStar',@(x) isa(x,'source'));
            inputs.addOptional('mmseStar',[],@(x) isa(x,'source'));
            inputs.addParamValue('telescope',[],@(x) isa(x,'telescopeAbstract'));
            inputs.addParamValue('pupil',true(sampling),@islogical);
            inputs.addParamValue('unit',[],@isnumeric);
            inputs.addParamValue('model','zonal',@ischar);
            inputs.addParamValue('zernikeMode',[],@isnumeric);
            inputs.addParamValue('noiseCovariance',[],@isnumeric);
            
            inputs.parse(sampling,diameter,atmModel,guideStar,varargin{:});
            
            obj.sampling   = inputs.Results.sampling;
            obj.diameter   = inputs.Results.diameter;
            obj.atmModel   = inputs.Results.atmModel;
            obj.guideStar  = inputs.Results.guideStar;
            obj.nGuideStar = length(obj.guideStar);
            obj.p_mmseStar = inputs.Results.mmseStar;    
            obj.nmmseStar  = length(obj.p_mmseStar);
            obj.pupil      = inputs.Results.pupil;   
            obj.unit       = inputs.Results.unit;   
            obj.model      = inputs.Results.model;
            obj.zernikeMode= inputs.Results.zernikeMode;
            obj.noiseCovariance   = inputs.Results.noiseCovariance;
            
            obj.guideStarListener = addlistener(obj,'guideStar','PostSet',@obj.resetGuideStar);
            obj.atmModel.wavelength = obj.p_mmseStar(1).wavelength;
            if isempty(inputs.Results.telescope)
                obj.tel = telescope(obj.diameter);
            else
                obj.tel = inputs.Results.telescope;
            end
            obj.log = logBook.checkIn(obj);
            
%             add(obj.log,obj,'Computing the covariance matrices')
            
            poolWasAlreadyOpen = true;
            
            switch obj.model
                
                case 'zonal'
                    
                    if matlabpool('size')==0
                        matlabpool('open')
                        poolWasAlreadyOpen = false;
                    end
                    
                    [obj.Cxx,obj.Cox] = ...
                        phaseStats.spatioAngularCovarianceMatrix(...
                        obj.sampling,obj.diameter,...
                        obj.atmModel,obj.guideStar,obj.p_mmseStar,'mask',obj.pupil);
                    obj.Coo = phaseStats.spatioAngularCovarianceMatrix(...
                        obj.sampling,obj.diameter,...
                        obj.atmModel,obj.p_mmseStar(1),'mask',obj.pupil);
            
                    if ~poolWasAlreadyOpen
                        matlabpool('close')
                    end
            
                case 'modal'
                    
                      obj.zernP    = zernike(obj.zernikeMode,obj.diameter,...
                          'resolution',obj.sampling);
                      obj.Cxx = zernikeStats.angularCovarianceAlt(obj.zernP,...
                          obj.atmModel,obj.guideStar,obj.guideStar);            
                      obj.Cox = { zernikeStats.angularCovarianceAlt(obj.zernP,...
                          obj.atmModel,obj.p_mmseStar,obj.guideStar) }; 
%                       obj.Cxx = cell2mat(zernikeStats.angularCovariance(obj.zernP,...
%                           obj.atmModel,obj.guideStar));            
%                       obj.Cox = zernikeStats.angularCovariance(obj.zernP,...
%                           obj.atmModel,obj.guideStar,obj.p_mmseStar)'; 
                      obj.Coo = zernikeStats.covariance(obj.zernP,obj.atmModel);
                    
                otherwise
                    
                    error('oomao:linearMMSE:linearMMSE',...
                        'The model is either zonal or modal')
                    
            end
            
%             add(obj.log,obj,'Computing the mmse builder')
            
            solveMmse(obj);
            
        end
        
        
        %% Destructor
        function delete(obj)
            if ~isempty(obj.log)
                checkOut(obj.log,obj)
            end
        end
        
        %% Set/Get mmseStar
        function set.mmseStar(obj,val)
            obj.p_mmseStar = val;
            add(obj.log,obj,'Computing the mmse/guide stars covariance matrix')
            obj.Cox = phaseStats.spatioAngularCovarianceMatrix(...
                obj.sampling,obj.diameter,...
                obj.atmModel,obj.guideStar,obj.p_mmseStar,'mask',obj.pupil);
            obj.Coo = phaseStats.spatioAngularCovarianceMatrix(...
                obj.sampling,obj.diameter,...
                obj.atmModel,obj.p_mmseStar(1),'mask',obj.pupil);
        end
        function val = get.mmseStar(obj)
            val = obj.p_mmseStar;
        end
        
        %% Set/Get noiseCovariance
        function set.noiseCovariance(obj,val)
            if ~isempty(val)
                if isempty(obj.zernNoiseCovariance)
                    zernTmp = zernike(1:66,obj.diameter);
                    obj.zernNoiseCovariance = zernTmp.noiseCovariance(obj.zernikeMode,obj.zernikeMode);
                end
                obj.nGuideStar = length(val);
                obj.p_noiseCovariance = repmat({obj.zernNoiseCovariance},1,obj.nGuideStar);
                obj.p_noiseCovariance = blkdiag( obj.p_noiseCovariance{:} );
                nMode = length(obj.zernikeMode);
                val = val(:);
                val = repmat( val , 1,  nMode )';
                obj.p_noiseCovariance = diag(val(:));%.*obj.p_noiseCovariance;
                solveMmse(obj);
                obj.Bmse = [];
            end
        end
        function val = get.noiseCovariance(obj)
            val = obj.p_noiseCovariance;
        end
        
        function lmmse(obj,fieldStar)
            %% LMMSE Linear Minimum Mean Square Error
            
            %             add(obj.log,obj,'Bayesian Minimum Mean Square Error computation in progress...')
            
            if nargin==1
                fun = @(x,y) obj.Coo - y*x';
                obj.Bmse = cellfun( fun , obj.Cox , obj.mmseBuilder , 'uniformOutput' , false );
            else
                Co1x = ...
                    phaseStats.spatioAngularCovarianceMatrix(...
                    obj.sampling,obj.diameter,...
                    obj.atmModel,obj.guideStar,fieldStar,'mask',obj.pupil);
                fun = @(x) obj.Coo + ...
                    obj.Cox{1}/obj.Cxx*obj.Cox{1}' - ...
                    x/obj.Cxx*obj.Cox{1}' - obj.Cox{1}/obj.Cxx*x';
                obj.Bmse = cellfun( fun , Co1x , 'uniformOutput' , false );
            end
            
        end
        
        function out = otf(obj,zRho)
            %% OTF Optical Transfer Function
                        
            if isempty(obj.Bmse)
                lmmse(obj)
            end
                        
            out   = zeros(size(zRho));
            rho   = abs(zRho);
            index = rho<obj.diameter;
            rho   = rho(index);
            
            rhoX = 2*real(zRho(index))/obj.diameter;
            rhoY = 2*imag(zRho(index))/obj.diameter;
            
            a = obj.Bmse{1};
            
            sf = a(1,1)*rhoX.^2 + a(2,2)*rhoY.^2 + (a(1,2)+a(2,1)).*rhoX.*rhoY;
            sf = 4*sf;
            
%             r0 = obj.atmModel.r0;
%             f0 = 1./obj.atmModel.L0;
%             fun = @(fx,fy,fh_rho) integrand(fx,fy,fh_rho,r0,f0);            
%             fc = 0.5*50/25.4;
%             sfFit = zeros(size(rho));
%             fxy = linspace(-fc,fc,1001);
%             [fx,fy] = meshgrid(fxy);
%             parfor k=1:numel(rho)
% %                sfFit(k) = dblquad( @(fx,fy) fun(fx,fy,rho(k)), -fc, fc , -fc , fc); 
%                sfFit(k) = trapz(fxy, trapz(fxy, fun(fx,fy,rho(k) ) ) ); 
%             end
%             sf = sf + phaseStats.structureFunction(rho,obj.atmModel) - sfFit;
% %             fun = @(f,x) phaseStats.spectrum(f,obj.atmModel).*(1 - besselj(0,2*pi*x.*f));
% %             sfFit = 4*pi*arrayfun( @(x) quadgk(@(f) fun(f,x),fc,Inf), rho(index));

            out(index) = otf(obj.tel,zRho(index)).*exp(-0.5*sf);
            
%             function out1 = integrand(fx,fy,rho_,r0,f0)
%                 
%                 f = hypot(fx,fy);
%                 out1 = (24.*gamma(6./5)./5).^(5./6).*...
%                     (gamma(11./6).^2./(2.*pi.^(11./3))).*...
%                     r0.^(-5./3);
%                 out1 = out1.*(f.^2 + f0.^2).^(-11./6).*(1 - besselj(0,2*pi*rho_.*f));
%                 
%             end
            
        end
        
        function out = image(obj,resolution,pixelScaleInMas)
            %% IMAGE 2D Point Spread Function
            %
            % psf = image(obj,resolution,pixelScaleInMas) computes the 2D
            % residual psf sampled with resolutionXresolution pixels of
            % size pixelScaleInMas; the psf is scaled such as its maximum
            % corresponds to the Strehl ratio
            
%             pxScaleAtNyquist = 0.25*obj.p_mmseStar.wavelength/obj.diameter;
%             imgLens = lens;
%             imgLens.nyquistSampling = ...
%                 pxScaleAtNyquist/(pixelScaleInMas*1e-3*constants.arcsec2radian);
%             n = resolution;
%             u = linspace(-1,1,n)*obj.diameter;
%             [x,y] = meshgrid(u);
%             z = x + 1i*y;
%             thisOtf = otf(obj,z);
%             src = source.*thisOtf*imgLens;
%             out = src.amplitude*obj.strehlRatio/max(src.amplitude(:));

              pixelScale = pixelScaleInMas*1e-3*constants.arcsec2radian/obj.p_mmseStar.wavelength;
              
              [fx,fy] = freqspace(resolution,'meshgrid');
              fx = pixelScale*fx*resolution/2;
              fy = pixelScale*fy*resolution/2;
            
              fc = 1;
              psd = fourierAdaptiveOptics.fittingPSD(fx,fy,fc,obj.atmModel);
              sf  = fft2(fftshift(psd))*pixelScale^2;
              sf  = 2*fftshift( sf(1) - sf );
            
              [rhoX,rhoY] = freqspace(resolution,'meshgrid');
              rhoX = 0.5*rhoX/pixelScale;
              rhoY = 0.5*rhoY/pixelScale;
%               rho  = hypot(rhoX,rhoY);

              [u,v] = freqspace(resolution,'meshgrid');
              fftPhasor = exp(1i.*pi.*(u+v)*0.5);

              thisOtf = otf(obj,rhoX+1i.*rhoY).*exp(-0.5*sf);
%               figure
%               mesh(rhoX,rhoY,thisOtf)
              out = real(ifftshift(ifft2(ifftshift(fftPhasor.*thisOtf))))/pixelScale^2;
              out = out./obj.tel.area;
%               out = out.*obj.strehlRatio./max(out(:));

        end
        
        function imagesc(obj,resolution,pixelScaleInMas)
            %% IMAGESC
            
            psf = image(obj,resolution,pixelScaleInMas);
            imagesc(psf)
            axis equal tight
            colorbar
        end
        
        function out = enSquaredEnergy(obj,eHalfSize)
            %% ENSQUAREDENERGY
            %
            % ee = enSquaredEnergy(obj,eHalfSize)
            
            a = 2*eHalfSize;
            out = quad2d(...
                @(o,r) r.*otf(obj,r.*exp(1i*o)).*...
                (sin(pi.*r.*cos(o).*a)./(pi.*r.*cos(o).*a)).*...
                (sin(pi.*r.*sin(o).*a)./(pi.*r.*sin(o).*a)), ...
                0,2*pi,0,obj.diameter).*a.*a;
        end
        
        function out = strehlRatio(obj)
            %% STREHLRATIO
            
            out = quad2d(...
                @(o,r) r.*otf(obj,r.*exp(1i*o)),0,2*pi,0,obj.diameter)./(pi*obj.diameter^2/4);
        end

        function map = get.varMap(obj)
            %% VARMAP Pupil error variance map
            
            if isempty(obj.Bmse)
                lmmse(obj)
            end
%             add(obj.log,obj,'Pupil error variance map computation in progress...')
            out = cellfun( @diag , obj.Bmse' , 'uniformOutput', false );
            out = cell2mat(out);
            map = zeros(obj.sampling,obj.sampling*obj.nmmseStar);
            mask = repmat( obj.pupil, 1 , obj.nmmseStar);
            map(mask) = out;
                
        end
       
        function out = get.var(obj)
            %% VAR Pupil error variance
            
            if isempty(obj.Bmse)
                lmmse(obj)
            end
%             add(obj.log,obj,'Pupil error variance computation in progress...')
            if iscell(obj.Bmse)
                out = cellfun( @trace , obj.Bmse , 'uniformOutput', true );
            else
                out = trace(obj.Bmse);
            end
            if strcmp(obj.model,'zonal')
                out = out/sum(obj.pupil(:));
            end
        end
        
        function map = get.rmsMap(obj)
            %% RMSMAP Pupil error rms map
            
            map = opd(obj,obj.varMap);
        end
        
        function out = get.rms(obj)
            %% RMS Pupil error rms
            
            out = opd(obj,obj.var);
        end        
        
        function out = get.rmsMas(obj)
            %% RMS Pupil error rms
            
            out = 1e3*constants.radian2arcsec*4*(sqrt(obj.var)/obj.p_mmseStar.waveNumber)/obj.diameter;
        end        

    end
       
    
    methods (Access=private)
        
        function out = opd(obj,val)
            %% OPD Phase to opd conversion
            
            out = sqrt(val);
            if ~isempty(obj.unit)
%                 add(obj.log,obj,sprintf('Rms in 1E%d meter',obj.unit))
                out = 10^-obj.unit*...
                    out*obj.atmModel.wavelength/(2*pi);
            end
        end
        
        function resetGuideStar(obj,varargin) % obj, src, event
            %% RESETGUIDESTAR
            
%             fprintf(' Resetting Cxx\n')
            switch obj.model
                
                case 'zonal'
                    
                    [obj.Cxx,obj.Cox] = ...
                        phaseStats.spatioAngularCovarianceMatrix(...
                        obj.sampling,obj.diameter,...
                        obj.atmModel,obj.guideStar,obj.p_mmseStar,'mask',obj.pupil);
                    
                case 'modal'
                    
                    obj.Cxx = zernikeStats.angularCovarianceAlt(obj.zernP,...
                        obj.atmModel,obj.guideStar,obj.guideStar);
                    obj.Cox = { zernikeStats.angularCovarianceAlt(obj.zernP,...
                        obj.atmModel,obj.p_mmseStar,obj.guideStar) };
                    
            end
            
            solveMmse(obj);
            obj.Bmse = [];
            
        end
        
        function solveMmse(obj)
            %% SOLVEMMSE
            
            m_mmseBuilder = cell(obj.nmmseStar,1);
            m_Cox = obj.Cox;
            m_Cxx = obj.Cxx;
            m_noiseCovariance = obj.p_noiseCovariance;
            parfor k=1:obj.nmmseStar
                m_mmseBuilder{k} = m_Cox{k}/(m_Cxx+m_noiseCovariance);
            end
            obj.mmseBuilder = m_mmseBuilder;

        end
        
    end
    
    methods (Static)
        
        function varargin = demoModal
            %% DEMOZONAL
            
            D = 10;
            n = 21;
            band = photometry.V;
            atm = atmosphere(photometry.V,0.15,60,'altitude',10e3);
            atm.wavelength = band;
%             gs = source('asterism',{[1,arcsec(60),0]},'wavelength',band);
            ss  = source('wavelength',band);
            pause(5)
            
%             mmse = linearMMSE(32,10,atm,gs,ss,...
%                 'model','modal','zernikeMode',2:3,...
%                 'unit',-9);
            zern = zernike(2:3,D,'resolution',n);
            %ZP = zern.p(zern.pupilLogical,:);
            
            gsRadius = 5:5:60;
            mmseTTRms = zeros(length(gsRadius),2);
            aniso   = zeros(1,length(gsRadius));
            m_log = logBook.checkIn();
            m_log.verbose = false;
            for k=1:length(gsRadius)
                fprintf(' --> gsRadius = %d arcsec\n', gsRadius(k) )
                gs = source('asterism',{[3,arcsec(gsRadius(k)),0]},'wavelength',band);
                mmseTT = linearMMSE(n,D,atm,gs,ss,...
                    'model','modal','zernikeMode',2:3,...
                    'unit',-9);
                mmseZonal = linearMMSE(n,D,atm,gs,ss,...
                    'pupil',zern.pupilLogical,'unit',-9);
                lmmse(mmseZonal)
                mmseTTRms(k,1) = mmseTT.rms;
                mmseTTRms(k,2) = mmseZonal.rms;%1e9*sqrt(trace(ZP'*mmseZonal.Bmse{1}*ZP))./ss.waveNumber;
                aniso(k) = zernikeStats.anisokinetism(zern,atm,gs(1),-9);
            end
            m_log.verbose = true;e
            
            figure
            plot(gsRadius,mmseTTRms,'.-',gsRadius,aniso,'.-')
            xlabel('GS radius [arcsec]')
            ylabel('WFE [nm]')
            legend('TT Stars','NGS Stars','TT anisok.')
            
            varargin = {atm, zern, mmseTT, mmseZonal};
            
        end
        
    end
    
end