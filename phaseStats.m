classdef phaseStats
    % PHASESTATS Phase statistics static class
    
    methods (Static)
        
        function  out = variance(atm)
            % VARIANCE Phase variance
            %
            % out = phaseStats.variance(atm) computes the phase variance
            % from an atmosphere object
            %
            % See also atmosphere
            
            L0r0ratio= (atm.L0./atm.r0).^(5./3);
            out   = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).*gamma(5./6)./(2.*pi.^(8./3))).*L0r0ratio;
            out = sum([atm.layer.fractionnalR0]).*out;
        end
        
        function out = covariance(rho,atm)
            % COVARIANCE Phase covariance
            %
            % out = phaseStats.covariance(rho,atm) computes the phase covariance from
            % the baseline rho and an atmosphere object
            %
            % See also atmosphere
            
            L0r0ratio= (atm.L0./atm.r0).^(5./3);
            cst      = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6)./(2.^(5./6).*pi.^(8./3))).*...
                L0r0ratio;
            out   = ones(size(rho)).*(24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).*gamma(5./6)./(2.*pi.^(8./3))).*L0r0ratio;
            index         = rho~=0;
            u             = 2.*pi.*rho(index)./atm.L0;
            out(index) = cst.*u.^(5./6).*besselk(5./6,u);
            out = sum([atm.layer.fractionnalR0]).*out;
        end
        
        function out = strctureFunction(rho,atm)
            % STRUCTUREFUNCTION Phase structure function
            %
            % out = phaseStats.structureFunction(rho,atm) computes the
            % phase structure function from the baseline rho and an
            % atmosphere object
            %
            % See also atmosphere
            
            out = 2.*(variance(atm)-covariance(rho,atm));
        end
        
        function out = spectrum(f,atm)
            % SPECTRUM Phase power spectrum density
            %
            % out = phaseStats.spectrum(f,atm) computes the phase power
            % spectrum density from the spatial frequency f and an
            % atmosphere object
            %
            % See also atmosphere
            
            out = (24.*gamma(6./5)./5).^(5./6).*...
                (gamma(11./6).^2./(2.*pi.^(11./3))).*...
                atm.r0.^(-5./3);
            out = out.*(f.^2 + 1./atm.L0.^2).^(-11./6);
            out = sum([atm.layer.fractionnalR0]).*out;
        end
        
        function out = covarianceMatrix(varargin)
            % COVARIANCEMATRIX Phase covariance matrix
            %
            % out = phaseStats.covarianceMatrix(rho1,atm) Computes the phase
            % auto-covariance matrix from the vector rho1 and an atmosphere
            % object
            %
            % out = phaseStats.covarianceMatrix(rho1,rho2,atm) Computes the phase
            % cross-covariance matrix from the vectors rho1 and rho2 and an
            % atmosphere object
            % Examples:
            % covariance matrix on a 1 metre square grid sampled on 16
            % pixels :
            % [x,y] = meshgrid((0:15)*1/15);
            % g = phaseStats.covarianceMatrix(complex(x,y),atm);
            % imagesc(g), axis square, colorbar
            % covariance matrix on a 1 meter square grid sampled on 16
            % pixels with the same grid but displaced of 1 meter:
            % [x,y] = meshgrid((0:15)*1/15);
            % z = complex(x,y);
            % g = phaseStats.covarianceMatrix(z,z+1,atm);
            % imagesc(g), axis square, colorbar
            %
            % See also atmosphere
            
            error(nargchk(2,3,nargin))
            rho1 = varargin{1}(:);
            if nargin==2
                atm  = varargin{2};
                rho  = abs(bsxfun(@minus,rho1,rho1.'));
            else
                rho2 = varargin{2}(:);
                atm  = varargin{3};
                rho  = abs(bsxfun(@minus,rho1,rho2.'));
            end
            [nRho,mRho] = size(rho);
            blockSize = 5000;
            if max(nRho,mRho)>blockSize % Memory gentle
                l  = floor(nRho/blockSize);
                le = nRho - l*blockSize;
                p  = floor(mRho/blockSize);
                pe = mRho - p*blockSize;
                rho = mat2cell( rho, ...
                    [ones(1,l)*blockSize, le],...
                    [ones(1,p)*blockSize, pe]);
                out = cell2mat( ...
                    cellfun(@(x) phaseStats.covariance(x,atm), rho, 'UniformOutput', false) );
            else % Memory intensive
                L0r0ratio= (atm.L0./atm.r0).^(5./3);
                cst      = (24.*gamma(6./5)./5).^(5./6).*...
                    (gamma(11./6)./(2.^(5./6).*pi.^(8./3))).*...
                    L0r0ratio;
                out   = ones(size(rho)).*(24.*gamma(6./5)./5).^(5./6).*...
                    (gamma(11./6).*gamma(5./6)./(2.*pi.^(8./3))).*L0r0ratio;
                index         = rho~=0;
                u             = 2.*pi.*rho(index)./atm.L0;
                out(index) = cst.*u.^(5./6).*besselk(5./6,u);
                out = sum(atm.layer.fractionnalR0).*out;
            end
        end
        
        function out = zernikeVariance(zern,atm,tel)
            % VARIANCE Zernike coefficients variance
            %
            % out = variance(modes,atmosphere,telescope) computes the
            % variance of Zernike coefficients from the modes, the
            % atmosphere object and the telescope object
            %
            % out = variance(zernike,atmosphere,telescope) computes the
            % variance of Zernike coefficients from the Zernike polynomials
            % object, the atmosphere object and the telescope object
            %
            % out = variance(modes,atmTel) computes the variance of Zernike
            % coefficients from the modes and the atmosphere&telescope
            % object
            %
            % out = variance(zernike,atmTel) computes the variance of
            % Zernike coefficients from the Zernike polynomials object and
            % the atmosphere&telescope object
            %
            % Example:
            % atm = atmosphere(photometry.V,0.15,30);
            % tel = telescope(10);
            % modes = 1:15;
            % figure
            % semilogy(modes,phaseStats.zernikeVariance(modes,atm,tel),'--.')
            % xlabel('Zernike modes')
            % ylabel('Variance [rd^2]')
            %
            % See also zernikepolynomials, atmosphere, telescope, AT
            
            
            if ~isa(zern,'zernike')
                zern = zernike(zern);
            end
            if nargin<3
                r0 = atm.r0;
                L0 = atm.L0;
                D  = atm.D;
            else
                r0 = atm.r0;
                L0 = atm.L0;
                D  = tel.D;
            end
            jv = zern.j;
            nv = zern.n;
            mv = zern.m;
            nv0   = nv;
            index = diff(nv)~=0;
            jv = [jv(index) jv(end)];
            mv = [mv(index) mv(end)];
            nv = [nv(index) nv(end)];
            nf  = length(nv);
            out = zeros(length(zern.j),1);
            
            for cpt = 1:nf
                
                j = jv(cpt);
                n = nv(cpt);
                m = mv(cpt);
                
                out(nv0==n,1) = zernCovCoef(r0,L0,D,j,j,n,m,n,m);
                
            end
            function out = zernCovCoef(r0,L0,D,i,j,ni,mi,nj,mj)
                if (mi==mj) && (rem(abs(i-j),2)==0 || ((mi==0) && (mj==0)))
                    if L0==Inf
                        if i==1 && j==1
                            out = Inf;
                        else
                            out = (gamma(11./6).^2.*gamma(14./3)./(2.^(8./3).*pi)).*(24.*gamma(6./5)./5).^(5./6).*...
                                (D./r0).^(5./3).*sqrt((ni+1).*(nj+1)).*(-1).^((ni+nj-mi-mj)./2).*...
                                newGamma(-5./6+(ni+nj)./2,...
                                [23./6+(ni+nj)./2 17./6+(ni-nj)./2 17./6+(nj-ni)./2]);
                        end
                    else
                        out = (4.*gamma(11./6).^2./pi.^(14./3)).*(24.*gamma(6./5)./5).^(5./6).*...
                            (L0./r0).^(5./3).*(L0./D).^2.*...
                            sqrt((ni+1).*(nj+1)).*(-1).^((ni+nj-mi-mj)./2).*...
                            UnParamEx4q2(0,ni+1,nj+1,11./6,pi.*D./L0);
                    end
                else
                    out = 0;
                end
                function out = newGamma(a,b)
                    % NEWGAMMA Computes the function defined by Eq.(1.18) in R.J. Sasiela's book :
                    % Electromagnetic Wave Propagation in Turbulence, Springer-Verlag.
                    % out = newGamma(a,b)
                    
                    out = prod(gamma(a))./prod(gamma(b));
                end
                function out = UnParamEx4q2(mu,alpha,beta,p,a)
                    % UNPARAMEX4Q2 Computes the integral given by the Eq.(2.33) of the thesis
                    % of R. Conan (Modelisation des effets de l'echelle externe de coherence
                    % spatiale du front d'onde pour l'Observation a Haute Resolution Angulaire
                    % en Astronomie, University of Nice-Sophia Antipolis, October 2000)
                    % http://www-astro.unice.fr/GSM/Bibliography.html#thesis
                    
                    a1 = [(alpha+beta+1)./2 (2+mu+alpha+beta)./2 (mu+alpha+beta)./2];
                    b1 = [1+alpha+beta 1+alpha 1+beta];
                    a2 = [(1-mu)./2+p 1+p p];
                    b2 = [1+(alpha+beta-mu)./2+p 1+(alpha-beta-mu)./2+p 1+(beta-alpha-mu)./2+p];
                    
                    out = (1./(2.*sqrt(pi).*gamma(p))).*(...
                        newGamma([a1 p-(mu+alpha+beta)./2],b1).*a.^(mu+alpha+beta).*...
                        pochammerSeries(3,5,a1,[1-p+(mu+alpha+beta)./2 b1 1],a.^2) + ...
                        newGamma([(mu+alpha+beta)./2-p a2],b2).*a.^(2.*p).*...
                        pochammerSeries(3,5,a2,[1-(mu+alpha+beta)./2+p b2 1],a.^2));
                    function out = pochammerSeries(p,q,a,b,z,tol,nmax)
                        % POCHAMMERSERIES Computes power series in Pochammer notation
                        % pochammerSeries(p,q,a,b,z)
                        % pochammerSeries(p,q,a,b,z,tol)
                        % pochammerSeries(p,q,a,b,z,[],nmax)
                        % pochammerSeries(p,q,a,b,z,tol,nmax)
                        
                        if (p==(q+1) && abs(z)<1) || (abs(z)==1 && real(sum(a)-sum(b))<0) || p<(q+1)
                            
                            if p==length(a) && q==length(b)
                                
                                switch nargin
                                    case 6
                                        nmax = 1e3;
                                    case 7
                                        if isempty(tol)
                                            tol = 1e-6;
                                        end
                                    otherwise
                                        tol = 1e-6;
                                        nmax = 1e3;
                                end
                                
                                out = zeros(size(z));
                                
                                indz = find(z==0);
                                if ~isempty(indz)
                                    out(indz) = 1;
                                end
                                
                                indnz = find(z~=0);
                                if ~isempty(indnz)
                                    z = z(indnz);
                                    ck = 1;
                                    step = Inf;
                                    k = 0;
                                    som = ck;
                                    while (k<=nmax) && (step>tol)
                                        ckp1 = prod(a+k).*z.*ck./prod(b+k);
                                        step = abs(abs(ck)-abs(ckp1));
                                        som = som + ckp1;
                                        k = k+1;
                                        ck = ckp1;
                                    end
                                    if step>tol
                                        warning('pochammerSeries','Maximum iteration reached before convergence')
                                    end
                                    out(indnz) = som;
                                end
                                
                            else
                                error('p and q must be the same length than vectors a and b, respectively')
                                
                            end
                            
                        else
                            error('This generalized hypergeometric function doesn''t converge')
                        end
                    end
                end
                
            end
        end
        
        function out = zernikeCovariance(zern,atm,tel)
            % COVARIANCE Zernike coefficients covariance
            %
            % out = phaseStats.covariance(modes,atmosphere,telescope)
            % computes the covariance matrix of Zernike coefficients from
            % the modes, the atmosphere object and the telescope object
            %
            % out = covariance(zernike,atmosphere,telescope) computes the
            % covariance matrix of Zernike coefficients from the Zernike
            % polynomials object, the atmosphere object and the telescope
            % object
            %
            % out = covariance(modes,atmTel) computes the covariance matrix
            % of Zernike coefficients from the modes and the
            % atmosphere&telescope object
            %
            % out = covariance(zernike,atmTel) computes the covariance
            % matrix of Zernike coefficients from the Zernike polynomials
            % object and the atmosphere&telescope object
            %
            % Example:
            % atm = atmosphere(photometry.V,0.15,30);
            % tel = telescope(10);
            % modes = 1:15;
            % figure
            % spy(phaseStats.zernikeCovariance(modes,atm,tel))
            %
            % See also zernikepolynomials, atmosphere, telescope
            
            
            if ~isa(zern,'zernike')
                zern = zernike(zern);
            end
            if nargin<3
                r0 = atm.r0;
                L0 = atm.L0;
                D  = atm.D;
            else
                r0 = atm.r0;
                L0 = atm.L0;
                D  = tel.D;
            end
            % -- covariance --
            n = zern.n;
            m = zern.m;
            [i,j]   = meshgrid(zern.j);
            [ni,nj] = meshgrid(n);
            [mi,mj] = meshgrid(m);
            index   = (mi==mj) & (rem(abs(i-j),2)==0 | ((mi==0) & (mj==0)));
            index   = logical(index - diag(diag(index)));
            i       = i(index);
            j       = j(index);
            ni      = ni(index);
            nj      = nj(index);
            mi      = mi(index);
            mj      = mj(index);
            
            nf       = sum(index(:));
            covValue = zeros(1,nf);
            for cpt = 1:nf
                covValue(cpt) = zernCovCoef(r0,L0,D,i(cpt),j(cpt),ni(cpt),mi(cpt),nj(cpt),mj(cpt));
            end
            out          = zeros(length(zern.j));
            out(index)   = covValue;
            % -- variance --
            jv = zern.j;
            nv = zern.n;
            mv = zern.m;
            nv0   = nv;
            index = diff(nv)~=0;
            jv = [jv(index) jv(end)];
            mv = [mv(index) mv(end)];
            nv = [nv(index) nv(end)];
            nf  = length(nv);
            variance = zeros(length(zern.j),1);
            
            for cpt = 1:nf
                
                j = jv(cpt);
                n = nv(cpt);
                m = mv(cpt);
                
                variance(nv0==n,1) = zernCovCoef(r0,L0,D,j,j,n,m,n,m);
                
            end
            % -- covariance + variance --
            out = out + diag(variance);
            function out = zernCovCoef(r0,L0,D,i,j,ni,mi,nj,mj)
                if (mi==mj) && (rem(abs(i-j),2)==0 || ((mi==0) && (mj==0)))
                    if L0==Inf
                        if i==1 && j==1
                            out = Inf;
                        else
                            out = (gamma(11./6).^2.*gamma(14./3)./(2.^(8./3).*pi)).*(24.*gamma(6./5)./5).^(5./6).*...
                                (D./r0).^(5./3).*sqrt((ni+1).*(nj+1)).*(-1).^((ni+nj-mi-mj)./2).*...
                                newGamma(-5./6+(ni+nj)./2,...
                                [23./6+(ni+nj)./2 17./6+(ni-nj)./2 17./6+(nj-ni)./2]);
                        end
                    else
                        out = (4.*gamma(11./6).^2./pi.^(14./3)).*(24.*gamma(6./5)./5).^(5./6).*...
                            (L0./r0).^(5./3).*(L0./D).^2.*...
                            sqrt((ni+1).*(nj+1)).*(-1).^((ni+nj-mi-mj)./2).*...
                            UnParamEx4q2(0,ni+1,nj+1,11./6,pi.*D./L0);
                    end
                else
                    out = 0;
                end
                function out = newGamma(a,b)
                    % NEWGAMMA Computes the function defined by Eq.(1.18) in R.J. Sasiela's book :
                    % Electromagnetic Wave Propagation in Turbulence, Springer-Verlag.
                    % out = newGamma(a,b)
                    
                    out = prod(gamma(a))./prod(gamma(b));
                end
                function out = UnParamEx4q2(mu,alpha,beta,p,a)
                    % UNPARAMEX4Q2 Computes the integral given by the Eq.(2.33) of the thesis
                    % of R. Conan (Modelisation des effets de l'echelle externe de coherence
                    % spatiale du front d'onde pour l'Observation a Haute Resolution Angulaire
                    % en Astronomie, University of Nice-Sophia Antipolis, October 2000)
                    % http://www-astro.unice.fr/GSM/Bibliography.html#thesis
                    
                    a1 = [(alpha+beta+1)./2 (2+mu+alpha+beta)./2 (mu+alpha+beta)./2];
                    b1 = [1+alpha+beta 1+alpha 1+beta];
                    a2 = [(1-mu)./2+p 1+p p];
                    b2 = [1+(alpha+beta-mu)./2+p 1+(alpha-beta-mu)./2+p 1+(beta-alpha-mu)./2+p];
                    
                    out = (1./(2.*sqrt(pi).*gamma(p))).*(...
                        newGamma([a1 p-(mu+alpha+beta)./2],b1).*a.^(mu+alpha+beta).*...
                        pochammerSeries(3,5,a1,[1-p+(mu+alpha+beta)./2 b1 1],a.^2) + ...
                        newGamma([(mu+alpha+beta)./2-p a2],b2).*a.^(2.*p).*...
                        pochammerSeries(3,5,a2,[1-(mu+alpha+beta)./2+p b2 1],a.^2));
                    function out = pochammerSeries(p,q,a,b,z,tol,nmax)
                        % POCHAMMERSERIES Computes power series in Pochammer notation
                        % pochammerSeries(p,q,a,b,z)
                        % pochammerSeries(p,q,a,b,z,tol)
                        % pochammerSeries(p,q,a,b,z,[],nmax)
                        % pochammerSeries(p,q,a,b,z,tol,nmax)
                        
                        if (p==(q+1) && abs(z)<1) || (abs(z)==1 && real(sum(a)-sum(b))<0) || p<(q+1)
                            
                            if p==length(a) && q==length(b)
                                
                                switch nargin
                                    case 6
                                        nmax = 1e3;
                                    case 7
                                        if isempty(tol)
                                            tol = 1e-6;
                                        end
                                    otherwise
                                        tol = 1e-6;
                                        nmax = 1e3;
                                end
                                
                                out = zeros(size(z));
                                
                                indz = find(z==0);
                                if ~isempty(indz)
                                    out(indz) = 1;
                                end
                                
                                indnz = find(z~=0);
                                if ~isempty(indnz)
                                    z = z(indnz);
                                    ck = 1;
                                    step = Inf;
                                    k = 0;
                                    som = ck;
                                    while (k<=nmax) && (step>tol)
                                        ckp1 = prod(a+k).*z.*ck./prod(b+k);
                                        step = abs(abs(ck)-abs(ckp1));
                                        som = som + ckp1;
                                        k = k+1;
                                        ck = ckp1;
                                    end
                                    if step>tol
                                        warning('pochammerSeries','Maximum iteration reached before convergence')
                                    end
                                    out(indnz) = som;
                                end
                                
                            else
                                error('p and q must be the same length than vectors a and b, respectively')
                                
                            end
                            
                        else
                            error('This generalized hypergeometric function doesn''t converge')
                        end
                    end
                end
                
            end
        end
        
        function out = zernikeResidualVariance(N,atm,tel)
            zern = zernike(1:N);
            if nargin<3
                r0 = atm.r0;
                L0 = atm.L0;
                D  = atm.D;
                aiVar = phaseStats.zernikeVariance(zern,atm);
            else
                r0 = atm.r0;
                L0 = atm.L0;
                D  = tel.D;
                aiVar = phaseStats.zernikeVariance(zern,atm,tel);
            end
            
            if isinf(L0)
                Delta1 = -(2.*gamma(11./6).^2./pi.^1.5).*(24.*gamma(6./5)./5).^(5./6).*...
                    (D./r0).^(5./3).*newGamma([-5./6,7./3],[23./6,17./6]);
                out = Delta1 - sum(aiVar(2:end));
            else
                out = phaseStats.variance(atm) - sum(aiVar);
            end
        end
        
        function aiaj = zernikeAngularCovariance(zern,atm,tel,src)
            % ZERNIKEANGULARCOVARIANCE Zernike coefficients angular
            % covariance
            %
            % aiaj = zernikeAngularCovariance(zern,atm,tel,src) computes
            % the covariance matrix between Zernike coefficients of Zernike
            % polynomials zern corresponding to wavefront propagating from
            % two sources src(1) and src(2) through the atmosphere atm and
            % to the telescope tel
            %
            % See also zernikepolynomials, atmosphere, telescope, source
            
            if src(1)==src(2)
                aiaj = phaseStats.zernikeCovariance(zern,atm,tel);
            else
                R   = tel.R;
                zs1 = src(1).height;
                zs2 = src(2).height;
                xSrc = tan(src(1).zenith).*cos(src(1).azimuth) - ...
                    tan(src(2).zenith).*cos(src(2).azimuth);
                ySrc = tan(src(1).zenith).*sin(src(1).azimuth) - ...
                    tan(src(2).zenith).*sin(src(2).azimuth);
                rhoSrcLayer = hypot(xSrc,ySrc);
                thetaSrcLayer = atan2(ySrc,xSrc);
                nMode = length(zern.j);
                znmj = mat2cell([zern.j;zern.n;zern.m],3,ones(1,nMode));
                znmj = repmat(znmj,nMode,1);
                znmi = znmj';
                index = triu(true(nMode));
                aiajFun = @ (znmi,znmj) ...
                    quadgk(@(x) integrand(x,znmi(1),znmi(2),znmi(3),znmj(1),znmj(2),znmj(3)), 0, Inf);
                aiaj = zeros(nMode);
                %                 tic
                aiaj(index) = cellfun(aiajFun,znmj(index),znmi(index));
                %                 toc
                aiaj = aiaj + triu(aiaj,1)';
                aiaj = bsxfun(@times,aiaj,(-1).^zern.m');
            end
            function out = integrand(x,zi,ni,mi,zj,nj,mj)
                krkr_mi = mi==0;
                krkr_mj = mj==0;
                out = 0;
                factor1 = sqrt((ni+1)*(nj+1)).*...
                    (-1).^(0.5.*(ni+nj)).*...
                    2.^(1-0.5.*(krkr_mi+krkr_mj));%.*...
                % (-1).^mj;
                for kLayer=1:atm.nLayer
                    a1l = R.*(1 - atm.layer(kLayer).altitude./zs1);
                    a2l = R.*(1 - atm.layer(kLayer).altitude./zs2);
                    denom = pi.*a1l.*a2l;
                    sl = atm.layer(kLayer).altitude.*rhoSrcLayer;
                    factor2 = (-1).^(1.5*(mi+mj)).*...
                        cos( ...
                        (mi+mj).*thetaSrcLayer + ...
                        pi.*( (1-krkr_mi).*((-1).^zi-1) + ...
                        (1-krkr_mj).*((-1).^zj-1) )./4 );
                    factor3 = (-1).^(1.5*abs(mi-mj)).*...
                        cos( ...
                        (mi-mj).*thetaSrcLayer + ...
                        pi.*( (1-krkr_mi).*((-1).^zi-1) - ...
                        (1-krkr_mj).*((-1).^zj-1) )./4 );
                    phasePSD = phaseStats.spectrum(0.5*x/pi,atm.slab(kLayer));
                    tripleBessel1 = besselj(mi+mj,sl.*x).*...
                        besselj(ni+1,a1l.*x).*...
                        besselj(nj+1,a2l.*x);
                    tripleBessel2 = besselj(abs(mi-mj),sl.*x).*...
                        besselj(ni+1,a1l.*x).*...
                        besselj(nj+1,a2l.*x);
                    out = out +  (factor1./denom).*(phasePSD./x).*...
                        ( factor2.*tripleBessel1 + factor3.*tripleBessel2 );
                end
                out = real(out);
            end
        end
        
        function out = zern_aiaj(zi,ni,mi,zj,nj,mj,atm,tel)
            out = 0;
            if (ni==0) && (nj==0)
                out = Inf;
            elseif rem(abs(zi-zj),2)==0 || (mi==0 && mj==0)
                red = pi*tel.D;
                out = 2/(pi*tel.R^2)*...
                    sqrt((ni+1)*(nj+1))*(-1)^((ni+nj-mi-mj)/2)*(mi==mj).*...
                    quadgk( @(x) phaseStats.spectrum(x,atm).*...
                    besselj(ni+1,red.*x).*besselj(nj+1,red.*x)./x,0,Inf);
            end
        end
        
        function out = zern_phiai(r,o,j,n,m,atm,tel)
            Inm = @(r) quadgk( @(f) phaseStats.spectrum(f,atm).*...
                besselj(n+1,pi.*f.*tel.D).*...
                besselj(m,2.*pi.*f.*r),0,Inf);
            out = arrayfun(Inm,r);
            krkr = m~=0;
            out = out.*2.*sqrt(n+1).*2.^(0.5*krkr).*(-1).^((n-m*krkr)/2).*...
                cos(m.*o+pi.*krkr.*((-1).^j-1)./4)./tel.R;
        end
        
        
        function map = zern_resid_var(N,atm,tel,r,o)
            zern = zernike(1:N);
            
            mapFun = @(rr,oo) quadgk(@(f)integrand(f,rr,oo),0,Inf);
            tic
            map = arrayfun( mapFun , r , o );
            toc
            function out = integrand(f,rr,oo)
                cum_i = 0;
                for i=1:N
                    cum_ip = 0;
                    for ip=1:N
                        flag = (zern.m(i)==zern.m(ip)) && ...
                            ( rem( abs( zern.j(i)-zern.j(ip) ),2 )==0 || ...
                            ( (zern.m(i)==0) && (zern.m(ip)==0) ) );
                        if flag
                            cum_ip = cum_ip + ...
                                sqrt(zern.n(ip)+1).*...
                                (-1).^((zern.n(ip)-zern.m(ip))/2).*...
                                %{
                                zernike.fun(zern.j(ip),zern.n(ip),zern.m(ip),rr,oo).*...
                                %}
                                besselj(zern.n(ip)+1,pi.*tel.D.*f);
                        end
                    end
                    cum_ip = cum_ip./(pi.*f.*tel.R);
                    deltaNot = 1 - zern.m(i)==0;
                    alpha = pi.*((-1).^zern.j(i)-1).*deltaNot./4;
%                     cum_ip = cum_ip - ...
%                         2.*(-1).^(zern.m(i).*(zern.m(i)==0)/2).*...
%                         2.^(deltaNot/2).*...
%                         cos(zern.m(i).*oo + alpha).*...
%                         besselj(zern.m(i),2.*pi.*f.*rr);
                    cum_i = cum_i + ...
                        sqrt(zern.n(i)+1).*...
                        (-1).^((zern.n(i)-zern.m(i))/2).*...
                        %{
                        zernike.fun(zern.j(i),zern.n(i),zern.m(i),rr,oo).*...
                        %}
                        besselj(zern.n(i)+1,pi.*tel.D.*f).*...
                        2.*cum_ip./tel.R;
                end
%                 out = phaseStats.spectrum(f,atm).*(2.*pi.*f + cum_i);
                out = phaseStats.spectrum(f,atm).*(cum_i);
            end
            
        end
        
    end
    
end