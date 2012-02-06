classdef fourierAdaptiveOptics < handle
    
    properties
        psf;
        paramListener;
        psfRootListener;
    end
    
    properties (SetObservable=true)
        tel;
        atm;
        src;
        noiseVariance;
        nActuator;
        psfResolution;
        psfPixelScaleInMas;
        loopGain;     % loop gain
        exposureTime; % exposure time
        latency;      % delay
        psfRoot = 1;
    end
    
    properties (Dependent)
        fc;
        fcInMas;
    end
    
    properties (Access=private)
        figHandle
    end
    
    methods 
        
        function obj = fourierAdaptiveOptics(tel,atm,nActuator,noiseVariance,g,T,tau)
            obj.tel           = tel;
            obj.atm           = atm;
            obj.noiseVariance = noiseVariance;
            obj.nActuator     = nActuator;
            if nargin>4
                obj.loopGain = g;
                obj.exposureTime = T;
                obj.latency = tau;
            end
            obj.paramListener = ...
                addlistener(obj,...
                {'noiseVariance','nActuator','tel',...
                'atm','psfResolution','psfPixelScaleInMas',...
                'loopGain','exposureTime','latency','src'},...
                'PostSet',@obj.resetPsf);
            obj.psfRootListener = addlistener(obj,'psfRoot','PostSet',@obj.resetPsfScale);
        end
        
        function out = get.fc(obj)
            out = 0.5*(obj.nActuator-1)/obj.tel.D;
        end
        function out = get.fcInMas(obj)
            out = obj.fc*obj.atm.wavelength*constants.radian2mas;
        end
            
        
        function out = fittingPSD(obj,fx,fy)
            %% FITTINGPSD Fitting error power spectrum density
            
            fc = obj.fc;
            out   = zeros(size(fx));
            index = abs(fx)>fc | abs(fy)>fc;
            f     = hypot(fx(index),fy(index));
            out(index) ...
                  = phaseStats.spectrum(f,obj.atm);
            out = out.*pistonFilter(obj,hypot(fx,fy));           
        end
        
        function out = noisePSD(obj,fx,fy)
            %% FITTINGPSD Fitting error power spectrum density
            
            fc = obj.fc;
            out   = zeros(size(fx));
            if obj.noiseVariance>0
                index = ~(abs(fx)>fc | abs(fy)>fc) & hypot(fx,fy)>0;
                f     = hypot(fx(index),fy(index));
                out(index) = obj.noiseVariance./...
                    ( 2*pi*f.*tools.sinc(0.5*fx(index)/fc).*tools.sinc(0.5*fy(index)/fc)).^2;
            out = out.*averageClosedLoopNoise(obj,fx,fy).*pistonFilter(obj,hypot(fx,fy));
            end
       end
        
        function out = aliasingPSD(obj,fx,fy)
            %% ALIASINGPSD Aliasing error power spectrum density
            
            fc    = obj.fc;
            out   = zeros(size(fx));
            index = ~(abs(fx)>fc | abs(fy)>fc);% & hypot(fx,fy)>1/obj.tel.D;
            pf = pistonFilter(obj,hypot(fx,fy));
            fx     = fx(index);
            fy     = fy(index);
%             figure,imagesc(ind),pause
            [fo,f] = cart2pol(fx,fy);
            n = 5;
            al = 0;
            nv = -n:n;
            for l=nv
                flx = fx - 2*l*fc;
                for m=nv
                    if (l~=0) && (m~=0)
                        fmy = fy - 2*m*fc;
                        flm = hypot(flx,fmy);
                        al = al + 0.25*sin(2*fo).^2.*(fx./fmy+fy./flx).^2.*...
                            phaseStats.spectrum(flm,obj.atm);
                    end
                end
            end
           l=0;
           flx = fx - 2*l*fc;
           ind = flx==0;
           nind = ~ind;
           nv(nv==0) = [];
           for m=nv
               fmy = fy - 2*m*fc;
               flm = hypot(flx,fmy);
               al(ind) = al(ind) + phaseStats.spectrum(flm(ind),obj.atm);
               al(nind) = al(nind) + 0.25*sin(2*fo(nind)).^2.*(fx(nind)./fmy(nind)+fy(nind)./flx(nind)).^2.*...
                   phaseStats.spectrum(flm(nind),obj.atm);
           end
           m = 0;
           fmy = fy - 2*m*fc;
           ind = fmy==0;
           nind = ~ind;
           for l=nv
               flx = fx - 2*l*fc;
               flm = hypot(flx,fmy);
               al(ind) = al(ind) + phaseStats.spectrum(flm(ind),obj.atm);
               al(nind) = al(nind) + 0.25*sin(2*fo(nind)).^2.*(fx(nind)./fmy(nind)+fy(nind)./flx(nind)).^2.*...
                   phaseStats.spectrum(flm(nind),obj.atm);
           end
           out(index) =  al.*averageClosedLoopAliasing(obj,fx,fy);
           out = out.*pf;
        end
        
        function out = servoLagPSD(obj,fx,fy)
            %% SERVOLAGPSD Servo-lag power spectrum density
            
            fc    = obj.fc;
            out   = zeros(size(fx));
            index = ~(abs(fx)>fc | abs(fy)>fc);
            pf = pistonFilter(obj,hypot(fx,fy));
            fx     = fx(index);
            fy     = fy(index);
            
            out(index) = phaseStats.spectrum(hypot(fx,fy),obj.atm).*averageClosedLoopRejection(obj,fx,fy);
            out = pf.*out;
        end
        
        function out = anisoplanatismPSD(obj,fx,fy)
            %% ANISOPLANATISM Anisoplanatism power spectrum density
            
            zLayer = [obj.atm.layer.altitude];
            fr0    = [obj.atm.layer.fractionnalR0];
            A = zeros(size(fx));
            for kLayer=1:obj.atm.nLayer
                red = 2*pi*zLayer(kLayer)*...
                    ( fx*obj.src.directionVector(1) + fy*obj.src.directionVector(2) );
                A  = A + fr0(kLayer)*( 1 - cos(red) );
            end
            out = pistonFilter(obj,hypot(fx,fy)).*A.*phaseStats.spectrum(hypot(fx,fy),obj.atm);
        end
        
        function out = powerSpectrumDensity(obj,fx,fy)
            %% POWERSPECTRUMDENSITY AO system power spectrum density
            
            out = fittingPSD(obj,fx,fy) + ...
                noisePSD(obj,fx,fy) + ...
                aliasingPSD(obj,fx,fy) + ...
                servoLagPSD(obj,fx,fy);
            if ~isempty(obj.src)
                out = out + anisoplanatismPSD(obj,fx,fy);
            end
        end
        
        function out = varFitting(obj)
            fc  = obj.fc;
            a = phaseStats.variance(obj.atm);
            b = dblquad( @(fx,fy) phaseStats.spectrum( hypot(fx,fy) , obj.atm ) , ...
                -fc,fc,-fc,fc);
            out = a - b;
        end
        
        function out = varServoLag(obj)
            fc  = obj.fc;
            out = quad2d( @(fx,fy) servoLagPSD(obj,fx,fy),-fc,fc,-fc,fc);
        end
        
        function out = varNoise(obj)
            fc  = obj.fc;
            out = quad2d( @(fx,fy) noisePSD(obj,fx,fy),-fc,fc,-fc,fc);
        end
        
        function out = varAo(obj,fLim)
            out = quad2d( @(fx,fy) powerSpectrumDensity(obj,fx,fy),-fLim,fLim,-fLim,fLim);
        end
        
        function varargout = image(obj,resolution,pixelScaleInMas)
            %% IMAGE Point spread function
            
            fprintf('Computing image plane ...\n')
            
            obj.psfResolution = resolution;
            obj.psfPixelScaleInMas = pixelScaleInMas;
            
            pixelScale = pixelScaleInMas*1e-3*constants.arcsec2radian/obj.atm.wavelength;

            [fx,fy] = freqspace(resolution,'meshgrid');
            fx = pixelScale*fx*resolution/2;
            fy = pixelScale*fy*resolution/2;
            
            psd = powerSpectrumDensity(obj,fx,fy);
            sf  = fft2(fftshift(psd))*pixelScale^2;
            sf  = 2*fftshift( sf(1) - sf );
            
            [rhoX,rhoY] = freqspace(resolution,'meshgrid');
            rhoX = 0.5*rhoX/pixelScale;
            rhoY = 0.5*rhoY/pixelScale;
            if all(abs(rhoX)<obj.tel.D)
                warning('oomao.fourierAdaptiveOptics.image','OTF is truncated: increase the resolution or decrease the pixel scale')
            end
%             rho  = hypot(rhoX,rhoY);
            telOtf = otf(obj.tel, rhoX+1i*rhoY);
            
            thisOtf = telOtf.*exp(-0.5*sf);
            
            [u,v] = freqspace(resolution,'meshgrid');
            fftPhasor = exp(1i.*pi.*(u+v)*0.5);
            
            obj.psf = real(ifftshift(ifft2(ifftshift(thisOtf.*fftPhasor))))/pixelScale^2;
            obj.psf = obj.psf/obj.tel.area;

            alpha = pixelScaleInMas*((0:resolution)-resolution/2);%(-resolution+1:2:resolution-1)/2;
            if isempty(obj.figHandle)
                obj.figHandle = figure;
            end
            figure(obj.figHandle)
            imagesc(alpha,alpha,obj.psf.^(1/obj.psfRoot))
            if any(abs(alpha)>obj.fcInMas)
                u = [-1 1 1 -1 -1]*obj.fcInMas;
                v = [-1 -1 1 1 -1]*obj.fcInMas;
                line(u,v,'color','r','linestyle',':')
            end
            axis xy equal tight
            xlabel('X axis [mas]')
            ylabel('Y axis [mas]')
            title(sprintf('Strehl: %4.2f%',sum(thisOtf(:))./sum(telOtf(:))*100))
            
            if nargout>0
                varargout{1} = obj.psf;
            end
%             pxScaleAtNyquist = 0.25/tel.D;
%             imgLens = lens;
%             imgLens.nyquistSampling = ...
%                 pxScaleAtNyquist/pixelScale;
            

        end
        
        function resetPsf(obj,varargin)
            if ~isempty(obj.psf)
                image(obj,obj.psfResolution,obj.psfPixelScaleInMas)
            end
        end
        
        function resetPsfScale(obj,varargin)
            if ~isempty(obj.psf)
                figure(obj.figHandle)
                h = findobj(obj.figHandle,'type','image');
                set(h,'Cdata',obj.psf.^(1/obj.psfRoot))
            end
        end
        
        function out = pistonFilter(obj,f)
            red = pi*obj.tel.D*f;
            out = 1 - 4*tools.sombrero(1,red).^2;
            
        end
        
        function out = closedLoopRejection(obj,nu)
            %% CLOSEDLOOPREJECTION Closed loop rejection transfer function
            
            out   = zeros(size(nu));
            index = nu~=0;
            nu   = nu(index);
            red = obj.loopGain.*tools.sinc(nu.*obj.exposureTime)./(2*pi*nu*obj.exposureTime);            
            out(index) = ....
                1./( 1 + red.^2 - 2.*red.*sin( 2*pi*nu*(obj.exposureTime+obj.latency) ) );
            
        end
        
        function E = averageClosedLoopRejection(obj,fx,fy)
            %% AVERAGECLOSEDLOOPREJECTION Atmosphere average closed loop rejection transfer function
            
            E = averageRejection(obj,fx,fy,@(nu)closedLoopRejection(obj,nu));
        end
        
        function out = closedLoopAliasing(obj,nu)
            %% CLOSEDLOOPALIASING Closed loop aliasing transfer function
            
            out   = ones(size(nu));
            index = nu~=0;
            nu   = nu(index);
            red = obj.loopGain.*tools.sinc(nu.*obj.exposureTime)./(2*pi*nu*obj.exposureTime);            
            out(index) = ....
                red.^2./( 1 + red.^2 - 2.*red.*sin( 2*pi*nu*(obj.exposureTime+obj.latency) ) );
            
        end
        
        function E = averageClosedLoopAliasing(obj,fx,fy)
            %% AVERAGECLOSEDLOOPALIASING Atmosphere average closed loop aliasing transfer function
            
            E = averageRejection(obj,fx,fy,@(nu)closedLoopAliasing(obj,nu));
        end
        
        function out = closedLoopNoise(obj,nu)
            %% CLOSEDLOOPALIASING Closed loop noise transfer function
            
            out   = ones(size(nu));
            index = nu~=0;
            nu   = nu(index);
            red = obj.loopGain.*tools.sinc(nu.*obj.exposureTime)./(2*pi*nu*obj.exposureTime);            
            out(index) = (red./tools.sinc(nu.*obj.exposureTime)).^2./...
                ( 1 + red.^2 - 2.*red.*sin( 2*pi*nu*(obj.exposureTime+obj.latency) ) );
            
        end
        
        function E = averageClosedLoopNoise(obj,fx,fy)
            %% AVERAGECLOSEDLOOPALIASING Atmosphere average closed loop noise transfer function
            
            E = averageRejection(obj,fx,fy,@(nu)closedLoopNoise(obj,nu));
        end
        
        function E = averageRejection(obj,fx,fy,fun)
            [vx,vy] = pol2cart([obj.atm.layer.windDirection],[obj.atm.layer.windSpeed]);
            fr0     = [obj.atm.layer.fractionnalR0];
            E = zeros(size(fx));
            for kLayer=1:obj.atm.nLayer
                nu = fx*vx(kLayer) + fy*vy(kLayer);
                E  = E + fr0(kLayer)*fun(nu);
            end
        end
        
    end
    
    methods (Static)
        
        function fao = gmtDemo
            wavelength = 1.6444e-6;
            d_atm = gmtAtmosphere(1);
            d_atm.wavelength = wavelength;
            d_atm.r0 = wavelength/(0.6/constants.radian2arcsec);
            resolution = 2^10;
            gmt = giantMagellanTelescope('resolution',resolution);
            fao = fourierAdaptiveOptics(gmt,d_atm,51,1.5,0.5,1/500,0);
            figure
            imagesc(image(gmt,resolution,1/wavelength/constants.radian2mas).^0.25)
            axis square
            image(fao,resolution,1)
        end
        
        function fao = psdDemo
            d_tel = telescope(8);
%             d_atm = atmosphere(photometry.V,15e-2,30,'windSpeed',10,'windDirection',0);
            d_atm = gmtAtmosphere(1);
            fao = fourierAdaptiveOptics(d_tel,d_atm,10,5,0.5,1e-3,1e-3);
            
            pixelScaleInMas = 2.5;
            pixelScale = pixelScaleInMas*1e-3*constants.arcsec2radian/fao.atm.wavelength;

            resolution = 128;
            [fx,fy] = freqspace(resolution,'meshgrid');
            fx = pixelScale*fx*resolution/2;
            fy = pixelScale*fy*resolution/2;
            f = hypot(fx,fy);
            
            alpha      = pixelScale*((0:resolution)-resolution/2);
            alpha      = alpha/fao.fc;

            figure
            subplot(2,2,1)
            imagesc(alpha,alpha,fittingPSD(fao,fx,fy))
            colorbar
            subplot(2,2,2)
            imagesc(alpha,alpha,noisePSD(fao,fx,fy))
            colorbar
            subplot(2,2,3)
            imagesc(alpha,alpha,aliasingPSD(fao,fx,fy))
            colorbar
            subplot(2,2,4)
            imagesc(alpha,alpha,servoLagPSD(fao,fx,fy))
            colorbar
       end
        
    end
    
end