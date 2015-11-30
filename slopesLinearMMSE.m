classdef slopesLinearMMSE < handle
    %% SLOPESLINEARMMSE MMSE wavefront estimation from slopes 
    %
    % mmse = linearMMSE(sampling,diameter,atmModel,guideStar)
    % mmse = linearMMSE(sampling,diameter,atmModel,guideStar,mmseStar)
    % mmse = linearMMSE(...,'pupil',pupilMask,'unit',unit)
    % mmse = linearMMSE(...,'model','modal','zernikeMode',zMode)
    % Results are given at the wavelength of the mmseStar
    %
    % Example:
    % sampling = 11;
    % diameter = 8;
    % atmModel = atmosphere(photometry.V,15e-2,60,'altitude',10e3);
    % gs = source('asterism',{[3,arcsec(30),0]},'wavelength',photometry.R);
    % ss = source('wavelength',photometry.H);
    % mmse = linearMMSE(sampling,diameter,atmModel,gs,ss,'pupil',utilities.piston(sampling,'type','logical'),'unit',-9);
    
    properties
        tag = 'slopesLinearMMSE';
        % pupil sampling
        sampling
        % pupil resolution
        resolution
        % pupil diameter
        diameter
        % slopes pupil mask
        slopesMask
        % wavefront pupil mask
        wavefrontMask
        % wavefront array size
        wavefrontSize;
        % atmosphere mode
        atmModel
        % guide star covariance matrix
        Cxx
        % mmse star covariance matrix
        Cox
        % guide star direction
        guideStar
        % Fourier spectrum resolution
        NF;
        % Fourier spectrum oversampling factor
        sf;
        % MINRES relative tolerance
        RTOL = 5e-2;
        % MINRES maximum number of iteration
        MAXIT = 200;
        % MINRES initial guess
        x0
        % warm start flag
        warmStart = false;
        % bilinear sparse operator
        B = [];
    end
    
    properties (Dependent)
        % mmse estimates directions
        mmseStar
        % noise variance
        noiseVar;
    end
    
    properties (Access=private)
        log
        p_mmseStar
        p_noiseVar = 0;
        fun;
        funp;
        c;
        wavefrontToMeter;
    end
 
    methods
        
        %% Constructor
        function obj = slopesLinearMMSE(wfs,tel,atmModel,guideStar,varargin)
            
            inputs = inputParser;
            inputs.addRequired('wfs',@(x) isa(x,'shackHartmann') );
            inputs.addRequired('tel',@(x) isa(x,'telescope'));
            inputs.addRequired('atmModel',@(x) isa(x,'atmosphere'));
            inputs.addRequired('guideStar',@(x) isa(x,'source'));
            inputs.addOptional('mmseStar',[],@(x) isa(x,'source'));
            inputs.addOptional('NF',1024, @isnumeric );
            inputs.addOptional('sf',4, @isnumeric );
           
            inputs.parse(wfs,tel,atmModel,guideStar,varargin{:});
            
            obj.resolution    = inputs.Results.wfs.lenslets.nLenslet;
            obj.diameter      = inputs.Results.tel.D;
            obj.sampling      = obj.diameter/obj.resolution;
            obj.wavefrontMask = ~inputs.Results.wfs.validActuator;
            obj.atmModel      = inputs.Results.atmModel;
            obj.guideStar     = inputs.Results.guideStar;
            obj.slopesMask    = repmat( repmat( inputs.Results.wfs.validLenslet(:) , 2, 1), 1, length(obj.guideStar) );
            obj.p_mmseStar    = inputs.Results.mmseStar;    
            obj.NF            = inputs.Results.NF;    
            obj.sf            = inputs.Results.sf;    
           
            obj.log = logBook.checkIn(obj);
            obj.log.verbose = false;
            add(obj.log,obj,'atmosphere wavelength set to mmse star wavelength!') 
            atmWavelength = obj.atmModel.wavelength;
            obj.atmModel.wavelength = obj.p_mmseStar(1).wavelength;
            
            if isinf(obj.guideStar(1).height)
                add(obj.log,obj,'Computing the covariance matrices for NGS')
                tic
                obj.Cxx = { slopestoSlopesCovariance(obj,0,0) };
                obj.Cox{1} = { { phaseToSlopesCovariance(obj,0,obj.atmModel,0,0,1) } };        
                toc
            else
                add(obj.log,obj,'Computing the covariance matrices for LGS')
                tic
                dv = [obj.guideStar.directionVector];
                ddv_x = bsxfun( @minus, dv(1,:), dv(1,:)' );
                ddv_y = bsxfun( @minus, dv(2,:), dv(2,:)' );
                m_Cxx = arrayfun( @(dx,dy) slopestoSlopesCovariance(obj, dx, dy), ...
                    ddv_x, ddv_y, 'UniformOutput', false );
                obj.Cxx = m_Cxx;
                
                obj.Cox{1} = cell( obj.atmModel.nLayer , 1 );
                for kLayer = 1:obj.atmModel.nLayer
%                     fprintf('. Layer #%d\n',kLayer);
                   atmLayer = slab(obj.atmModel,kLayer);
                    g   = 1 - atmLayer.layer.altitude./[obj.guideStar.height];
                    pad = ceil( 0.5*obj.resolution*(1-g)./g );
                    deltaSrc = bsxfun( @minus ,dv , obj.p_mmseStar.directionVector );
                    deltaSrc = atmLayer.layer.altitude.*deltaSrc;
                    m_Cox =...
                        arrayfun( @(x,dx,dy,gl) phaseToSlopesCovariance(obj,x,atmLayer,dx,dy,gl) ,...
                        pad , deltaSrc(1,:), deltaSrc(2,:), g, ...
                        'UniformOutput', false );
                    obj.Cox{1}{kLayer} = m_Cox;
                    NI = obj.resolution+1;
                    NP = NI + 2*pad(1);
                    obj.B{kLayer} = tools.bilinearSparseInterpolator(NI,NP,1,g(1));
                end
                obj.B = cell2mat( obj.B );
                toc
            end
            
            obj.fun = @(x) mtimes4squareBlocks(obj.Cxx,x,wfs.validLenslet(:));
            obj.funp = @(x) mtimes4precond(obj.Cxx,x,obj.slopesMask);
            obj.c   = zeros(obj.resolution^2*2*length(obj.guideStar),1 );
            obj.wavefrontToMeter = ...
                obj.guideStar(1).wavelength*obj.atmModel.wavelength/2/pi/(obj.sampling*wfs.lenslets.fftPad)/wfs.lenslets.nyquistSampling;
            obj.wavefrontSize    = ones(1,2)*(obj.resolution+1);
            obj.x0               = zeros(size(obj.c));
            obj.log.verbose = true;
            obj.atmModel.wavelength = atmWavelength;
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
            phaseToSlopesCovariance(obj);
        end
        function val = get.p_mmseStar(obj)
            val = obj.p_mmseStar;
        end
        
        %% Set/Get noiseVar
        function set.noiseVar(obj,val)
            obj.p_noiseVar = val;
            slopestoSlopesCovariance(obj);
        end
        function val = get.noiseVar(obj)
            val = obj.p_noiseVar;
        end
                
        %% Wavefront reconstruction
        function out = mtimes(obj,wfs)
            if isa(wfs,'shackHartmann')
                obj.c( obj.slopesMask ) = wfs.slopes(:);
            else
                obj.c( obj.slopesMask ) = wfs;
            end
            [yy,~] = minres(obj.fun,obj.c(:),obj.RTOL,obj.MAXIT,[],[],obj.x0*obj.warmStart);
            obj.x0 = yy;
            m_Cox = obj.Cox{1};
            n = length(m_Cox);
            m = length(m_Cox{1});
            out = cell( n , 1 );
            ns = obj.resolution^2*2;
%             nFig = 100 + randi(100);
            for kn=1:n
                out{kn} = zeros( m_Cox{kn}{1}{1}.nRow^2 , 1 );
                u  = 1:ns;
                for km=1:m
                    yyy = yy(u);
                    out{kn} = out{kn} + ...
                        m_Cox{kn}{km}{1}*yyy(1:end/2) + m_Cox{kn}{km}{2}*yyy(1+end/2:end);
                    u = u + ns;
                end
%                 figure(nFig)
%                 subplot(1,n,kn)
%                 imagesc( reshape( out{kn} , ones(1,2)*m_Cox{kn}{1}{1}.nRow ) )
%                 axis square 
%                 colorbar
            end
            out = cell2mat( out );
            if ~isempty(obj.B)
                out = obj.B*out;
            end
            if length(out)>prod(obj.wavefrontSize)
                out(obj.wavefrontMask) = [];
            else
                if numel(out)==numel(obj.wavefrontMask)
                    out(obj.wavefrontMask) = 0;
                end
                out = reshape(out,obj.wavefrontSize);
            end
            out = out*obj.wavefrontToMeter;
        end
        function out = estimation(obj,wfs)
            % [yy,flag,relres,iter,resvec] = minres(fun,c,1e-3,50);
            if isa(wfs,'shackHartmann')
                obj.c( obj.slopesMask ) = wfs.slopes;
            else
                obj.c( obj.slopesMask ) = wfs;
            end
            yy = minres(obj.fun,obj.c,obj.RTOL,obj.MAXIT,[],[],obj.x0);
            if obj.warmStart
                obj.x0 = yy;
            end
            out = obj.Cox{1}*yy(1:end/2) + obj.Cox{2}*yy(1+end/2:end);
            if length(out)>prod(obj.wavefrontSize)
                out(obj.wavefrontMask) = [];
            else
                out(obj.wavefrontMask) = 0;
                out = reshape(out,obj.wavefrontSize);
            end
            out = out*obj.wavefrontToMeter;
        end
    
        function varargout = imagesc(obj)
%             n = length(obj.Cxx);
            m_Cxx = cellfun( ...
                @(x) [x{1,1}.elements, x{1,2}.elements x{2,1}.elements, x{2,2}.elements] ,...
                obj.Cxx,'UniformOutput',false);
            m_Cxx = cell2mat( m_Cxx );
            m_Cox = cellfun( ...
                @(x) [x{1,1}.elements, x{1,2}.elements] ,...
                obj.Cox{1},'UniformOutput',false);
%             m_Cox = cell2mat( m_Cox );
%             subplot(2*n+1,2*n,[1,4*n^2])
%             imagesc(m_Cxx)
%             axis square
%             colorbar
%             subplot(2*n+1,2*n,[1+4*n^2,(2*n+1)*(2*n)])
%             imagesc(m_Cox)
%             axis equal tight
%             colorbar
            if nargout>0
                varargout{1} = m_Cxx;
                varargout{2} = m_Cox;
            end
        end
    end

    methods (Access=private)
        
        function  m_Cxx = slopestoSlopesCovariance(obj,deltaSrc_x,deltaSrc_y)
            %% slopes-to-slopes covariance matrix
%             fprintf(' ==> Computing slopes-to-slopes covariance matrix ...\n');
            nLenslet = obj.resolution;
            atm = obj.atmModel;
            lambda = atm.wavelength;
            d = obj.diameter/nLenslet;
            
            covxx = zeros( obj.NF);
            covyy = zeros( obj.NF);
            cov   = zeros( obj.NF);
            
            nm = ones(1,2)*nLenslet;
            b0 = obj.NF/2+1;
            b  = (1-nLenslet:nLenslet-1)*obj.sf + b0;
            
            [fx0,fy0] = freqspace(obj.NF,'meshgrid');
            
%             tic
%             fprintf(' . Layer #  ');
            for kLayer = 1:obj.atmModel.nLayer
                
%                 fprintf('\b\b%2d',kLayer);
                atmLayer = slab(obj.atmModel,kLayer);
                h = atmLayer.layer.altitude;
                g = 1-h/obj.guideStar(1).height;
                
                lf = obj.sf/(d*2*g);
                fx = lf*fx0;
                fy = lf*fy0;
            
                delta = 2*lf/obj.NF;
                phasor = exp( 2*1i*pi.*h*( deltaSrc_x*fx + deltaSrc_y*fy ) );
                spectrum = @(fx,fy,u,v) lambda.^2*(fx.*u(1) + fy.*u(2)).*(fx.*v(1) + fy.*v(2)).*...
                    delta.^2.*phaseStats.spectrum(hypot(fx,fy),atmLayer).*...
                    (tools.sinc(g*d*fx).*tools.sinc(g*d*fy)).^2.*phasor;
            
                covxx = covxx + fft2( fftshift( spectrum(fx,fy,[1,0],[1,0]) ) );
                covyy = covyy + fft2( fftshift( spectrum(fx,fy,[0,1],[0,1]) ) );
                cov   = cov   + real( fftshift( fft2( fftshift( spectrum(fx,fy,[0,1],[1,0]) ) ) ) );

            end
%             fprintf('\n');
            
            covxx(1) = covxx(1) + obj.p_noiseVar;
            covxx = real( fftshift(  covxx ) );
            T = toeplitzBlockToeplitz( nm, nm, covxx(b,b) );
            m_Cxx{1,1} = T;
            
            covyy(1) = covyy(1) + obj.p_noiseVar;
            covyy = real( fftshift( covyy ) );
            T = toeplitzBlockToeplitz( nm, nm, covyy(b,b) );
            m_Cxx{2,2} = T;
            
            T = toeplitzBlockToeplitz( nm, nm, cov(b,b) );
            m_Cxx{1,2} = T;
            m_Cxx{2,1} = T';
            
%             elapsedTime = toc;
%             fprintf(' ==> slopes-to-slopes covariance matrix computed in %5.2fs\n',elapsedTime);
        end
        
        function m_Cox = phaseToSlopesCovariance(obj,pad,m_atm,deltaSrc_x,deltaSrc_y, gl)
            %% phase-to-slopes covariance matrix
%             fprintf(' ==> Computing phase-to-slopes covariance matrix ...\n');
            alpha = 1;
            [fx,fy] = freqspace(obj.NF,'meshgrid');
            nLenslet = obj.resolution;
            lambda = m_atm.wavelength;
            d = obj.sampling;
            lf = obj.sf/(d*2*gl);
            fx = lf*alpha*fx;
            fy = lf*alpha*fy;
            delta = 2*lf*alpha/obj.NF;
            phasor = exp( 2*1i*pi.*( (deltaSrc_x+0.5*d)*fx + (deltaSrc_y+0.5*d)*fy ) );
            spectrum2 = @(u) lambda.*1i*(fx.*u(1) + fy.*u(2)).*...
                delta.^2.*phaseStats.spectrum(hypot(fx,fy),m_atm).*...
                tools.sinc(gl*d*fx).*tools.sinc(gl*d*fy).*phasor./gl;
            
            nm = [nLenslet+1+2*pad nLenslet];
            b0 = obj.NF/2+1;
            b = (1-nLenslet-pad:nLenslet+pad)*obj.sf + b0;
            
%             tic
            covx  = fftshift(real( fft2( fftshift( spectrum2([1,0]) ) ) ) );
            covy  = fftshift(real( fft2( fftshift( spectrum2([0,1]) ) ) ) );
            m_Cox{1} = toeplitzBlockToeplitz( nm, nm, covx(b,b) );
            m_Cox{2} = toeplitzBlockToeplitz( nm, nm, covy(b,b) );
%             elapsedTime = toc;
%             fprintf(' ==> phase-to-slopes covariance matrix computed in %5.2fs\n',elapsedTime);
        end
        
    end
    
end

function out = mtimes4squareBlocks( OO, w, mask )

n = length(OO);
nw = length(w);
u = 1:nw/n;
out = zeros( size(w) );
if nargin>2
    for ii=1:n
        v = 1:nw/n;
        for jj=1:n
            ww = w(v);
            x = ww(1:end/2);
            y = ww(1+end/2:end);
            O = OO{ii,jj};
            out(u) = out(u) + [...
                maskedMtimes(O{1,1},x,mask) + maskedMtimes(O{1,2},y,mask)
                maskedMtimes(O{2,1},x,mask) + maskedMtimes(O{2,2},y,mask)
                ];
            v = v + nw/n;
        end
        u = u + nw/n;
    end
else
    for ii=1:n
        v = 1:nw/n;
        for jj=1:n
            ww = w(v);
            x = ww(1:end/2);
            y = ww(1+end/2:end);
            O = OO{ii,jj};
            out(u) = out(u) + [...
                mtimes(O{1,1},x) + mtimes(O{1,2},y)
                mtimes(O{2,1},x) + mtimes(O{2,2},y)
                ];
            v = v + nw/n;
        end
        u = u + nw/n;
    end
end

end

function out = mtimes4precond( OO, w, mask)

x = w(1:end/2);
y = w(1+end/2:end);
O = OO{1};
A = O{1,1};
B = O{2,2};
out = [medfilt1( A\x ) ; medfilt1( B\y )];

end