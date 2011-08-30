classdef fourierAdaptiveOptics < handle
    
    properties
        psf;
        paramListener;
    end
    
    properties (SetObservable=true)
        tel;
        atm;
        noiseVariance;
        nActuator;
        psfResolution;
        psfPixelScaleInMas;
    end
    
    properties (Dependent)
        fc;
        fcInMas;
    end
    
    properties (Access=private)
        figHandle
    end
    
    methods 
        
        function obj = fourierAdaptiveOptics(tel,atm,nActuator,noiseVariance)
            obj.tel           = tel;
            obj.atm           = atm;
            obj.noiseVariance = noiseVariance;
            obj.nActuator     = nActuator;
            obj.paramListener = ...
                addlistener(obj,...
                {'noiseVariance','nActuator','tel',...
                'atm','psfResolution','psfPixelScaleInMas'},...
                'PostSet',@obj.resetPsf);
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
            end
            out = out.*pistonFilter(obj,hypot(fx,fy));
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
           out(index) =  al;
           out = out.*pf;
        end
        
        function varargout = image(obj,resolution,pixelScaleInMas)
            %% OMAGE Point spread function
            
            obj.psfResolution = resolution;
            obj.psfPixelScaleInMas = pixelScaleInMas;
            
            pixelScale = pixelScaleInMas*1e-3*constants.arcsec2radian/obj.atm.wavelength;

            [fx,fy] = freqspace(resolution,'meshgrid');
            fx = pixelScale*fx*resolution/2;
            fy = pixelScale*fy*resolution/2;
            
            psd = fittingPSD(obj,fx,fy) + noisePSD(obj,fx,fy) + aliasingPSD(obj,fx,fy);
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
            imagesc(alpha,alpha,obj.psf.^0.25)
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
        
        function out = pistonFilter(obj,f)
            red = pi*obj.tel.D*f;
            out = 1 - 4*tools.sombrero(1,red).^2;
            
        end
    end
    
    methods (Static)
        
        function fao = gmtDemo
            wavelength = 1.6444e-6;
            d_atm = gmtAtmosphere(1);
            d_atm.wavelength = wavelength;
            d_atm.r0 = wavelength/(0.6/constants.radian2arcsec);
            resolution = 88*9;
            gmt = giantMagellanTelescope('resolution',resolution);
            fao = fourierAdaptiveOptics(gmt,d_atm,51,1.5);
            image(fao,resolution,6/9)
        end
        
        function fao = psdDemo
            d_tel = telescope(8);
            d_atm = atmosphere(photometry.V,15e-2,30);
            fao = fourierAdaptiveOptics(d_tel,d_atm,10,5);
            
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
            imagesc(pistonFilter(fao,hypot(fx,fy)))
            colorbar
       end
        
    end
    
end