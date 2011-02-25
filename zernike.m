classdef zernike < telescopeAbstract
    %% Zernike polynomials
    %
    % obj = zernike(j) creates a Zernike polynomials object from the
    % modes vector. The radius is set to unity.
    %
    % obj = zernike(j,D) creates a Zernike polynomials object from the
    % modes vector and the telescope diameter D. 
    %
    % obj = zernike(...,'resolution',nPixel) creates a Zernike polynomials
    % object from the modes vector and the the number of pixel across the
    % polynomials diameter
    %
    % obj = zernike(...,'radius',r,'angle',o) creates a Zernike polynomials
    % object from the modes vector computed on the given radius and angle.
    % The radius must be normalized to one.
    %
    % obj = zernike(...,'resolution',nPixel,'pupil',pupilMask) creates a
    % Zernike polynomials object from the modes vector, the the number of
    % pixel across the polynomials diameter and the pupil mask
    %
    % Example:
    % zern = zernike(1:6,8,'resolution',32);
    % %Computes Zernike polynomials from piston to astigmastism on an 8m
    % pupil sampled on a 32x32 grid 
    % Projection of Zernike modes onto themselves
    % zern = zern.\zern.p;
    % disp(zern.c); % the coefficients
    % Least square fit of the zernike modes by themselves
    % zern = zern\zern.p;
    % disp(zern.c); % the coefficients
    
    properties
        % mode number
        j;
        % radial order
        n;
        % azimuthal frequency
        m;
        % normalized radius
        r;
        % angle
        o;
        % Noll normalisation
        nollNorm;
        % coefficients in meter
        c;
        % lexicographic ordering of frames (default: true)
        lex = true;
        %         % Projection matrix of a Zernike portion on the same basis defined
        %         % on the portion
        %         P
%         pupil;
        % zernike tag
        tag = 'ZERNIKE POLYNOMIALS';
    end
    
    properties (Dependent)%, SetAccess=private)
        % polynomials
        p;
        % # of modes
        nMode;
        % Zernike X derivative matrix
        xDerivativeMatrix
        % Zernike Y derivative matrix
        yDerivativeMatrix
        % Zernike X derivative
        xDerivative
        % Zernike Y derivative
        yDerivative
        % telescope pupil mask
        pupil;
        % alias to polynomials p
        modes
        % noise covariance matrix (Shack-Hartmann wavefront sensor)
        noiseCovariance 
    end
    
    properties (Access=private)
        p_p;
    end
    
    methods
        
        %% Constructor
        function obj = zernike(j,varargin)
            error(nargchk(1,16,nargin));
            p = inputParser;
            p.addRequired('j', @isnumeric);
            p.addOptional('D', 2, @isnumeric);
            p.addParamValue('obstructionRatio', 0, @isnumeric);
            p.addParamValue('resolution', [], @isnumeric);
            p.addParamValue('radius', [], @isnumeric);
            p.addParamValue('angle', [], @isnumeric);
            p.addParamValue('pupil', [], @(x) isnumeric(x) || islogical(x));
            p.addParamValue('fieldOfViewInArcmin', [], @isnumeric);
            p.addParamValue('unitNorm', false, @islogical);
            p.addParamValue('logging', true, @islogical);
            p.parse(j,varargin{:});
            obj = obj@telescopeAbstract(p.Results.D,...
                'obstructionRatio',p.Results.obstructionRatio,...
                'resolution',p.Results.resolution,...
                'fieldOfViewInArcmin',p.Results.fieldOfViewInArcmin);
%             obj.pupil = p.Results.pupil;
            obj.j = p.Results.j;
            obj.r = p.Results.radius;
            obj.o = p.Results.angle;
            [obj.n,obj.m] = findNM(obj);
            obj.nollNorm = (sqrt((2-(obj.m==0)).*(obj.n+1)))';
%             obj.nollNorm(1)=[];
            if isempty(obj.r)
                [obj.r,obj.o] = utilities.cartAndPol(p.Results.resolution,'output','polar');
            else
                obj.resolution = length(obj.r);
            end
            obj.p_p = polynomials(obj,p.Results.unitNorm);
            obj.c = ones(length(obj.j),1);
            obj.log = logBook.checkIn(obj);
            if obj.log.writable
                display(obj)
            end
        end
        
        %% Destructor
        function delete(obj)
            if ~isempty(obj.log)
                checkOut(obj.log,obj)
            end
        end

        %% Get and Set the pupil
        function pupil = get.pupil(obj)
            pupil = obj.p_pupil;
            if isempty(pupil) && ~isempty(obj.resolution)
                pupil = obj.r>=obj.obstructionRatio & obj.r<=1;
                obj.p_pupil = double(pupil);
            end
        end
        
        %% Get modes
        function modes = get.modes(obj)
            modes = obj.p_p;
        end
        
        function out = mtimes(obj,c)
            %% MTIMES Zernike multiplication
            %
            % v = obj*c multiply the Zernike modes in obj by the
            % vector c
            
            out = obj.modes*c;
        end

        
        function display(obj)
            %% DISPLAY Display object information
            %
            % display(obj) prints information about the atmosphere+telescope object
            
            fprintf('___ %s ___\n',obj.tag)
            fprintf(' . %d modes\n',obj.nMode)
            fprintf('----------------------------------------------------\n')
            
        end

        function obj = minus(obj,val)
            %% - Modes cancellation
            %
            % obj = obj - cancelMode removes the polynomials matching 
            % cancelMode from the Zernike basis
            
            index = ismember(obj.j,val);
            obj.j(index) = [];
            obj.n(index) = [];            
            obj.m(index) = [];            
            obj.c(index,:) = [];     
            if ~isempty(obj.p_p)
                obj.p_p(:,index) = [];
            end
        end

        function obj = plus(obj,val)
            %% - Modes expansion
            %
            % obj = obj + newMode adds polynomials matching 
            % newMode to the Zernike basis
            
            newObj = zernike(val,obj.D,...
                'resolution',obj.resolution,...
                'radius',obj.r,...
                'angle',obj.o,...
                'pupil',obj.pupil,...
                'logging',false);
            j = [obj.j newObj.j];
            n = [obj.n newObj.n];
            m = [obj.m newObj.m];
            p = [obj.p_p newObj.p_p];
            [obj.j,index] = sort(j);
            obj.n = n(index);
            obj.m = m(index);
            obj.p_p = p(:,index);
        end
        
        function obj = ldivide(obj,val)
            %% .\ Orthogonal projection onto the Zernike polynomials
            %
            % obj = obj.\phaseMap projects the 2D phase map onto
            % the polynomials of the zernike object and store the
            % projection coefficients into the object coefficients vector
            
            if isa(val,'source')
                nVal = length(val);
                phaseMap = zeros(length(val(1).phase(:)),nVal);
                for kVal=1:nVal
                    phaseMap(:,kVal) = val(kVal).phase(:)/val(kVal).waveNumber;
                end
            else
                phaseMap = val;
            end
            if size(obj.p_p,1)==size(phaseMap,1)
                obj.c = obj.p_p'*phaseMap;
            else
                obj.c = obj.p_p'*utilities.toggleFrame(phaseMap,2);
            end
            obj.c = bsxfun( @times , obj.c , 1./diag(obj.p_p'*obj.p_p) );
        end
        
        function obj = mldivide(obj,val) 
            %% \ Least square fit to the Zernike polynomials
            %
            % out = obj\phaseMap fits to the least square the 2D phase map
            % onto the polynomials of the zernike object and store the
            % projection coefficients into the object coefficients vector
            
            if isa(val,'shackHartmann')
                u = val.validLenslet;
                obj = zernike(1:obj.j(end),...
                    'resolution',val.lenslets.nLenslet,...
                    'pupil',double(u));
                val.zern2slopes = [obj.xDerivative(u,:);obj.yDerivative(u,:)];
                obj.c = val.zern2slopes\val.slopes;
                obj = obj - 1;
            else
                if isa(val,'source')
                    nVal = length(val);
                    phaseMap = zeros(length(val(1).phase(:)),nVal);
                    for kVal=1:nVal
                        phaseMap(:,kVal) = val(kVal).phase(:)/val(kVal).waveNumber;
                    end
                else
                    phaseMap = val;
                end
                if size(obj.p_p,1)==size(phaseMap,1)
                    obj.c = obj.p_p\phaseMap;
                else
                    obj.c = obj.p_p\utilities.toggleFrame(phaseMap,2);
                end
            end
        end
        
        function obj = times(obj,factor)
            %% * Zernike coefficients multiplication
            %
            % obj = obj*factor multiplies the Zernike coefficients by
            % factor; factor is either a scalar or a vector the same size
            % than Zernike coefficients
            
            obj.c = obj.c.*factor;
            
        end
        
        function varargout = uminus(obj) 
            %% uminus Coefficient negation
            %
            % -obj multiplies the Zernike coefficients by -1
            
            obj.c = -obj.c;
            if nargout>0
                varargout{1} = obj;
            end
        end
        
        %% get p
        function out = get.p(obj)
            out = obj.p_p;
            if ~obj.lex
                out = utilities.toggleFrame(out,3);
            end
        end
        
        %% nMode
        function out = get.nMode(obj)
            out = length(obj.j);
        end
        
        function relay(obj,src)
            %% RELAY zernike to source relay
            %
            % relay(obj,srcs) writes the zernike amplitude and phase into
            % the properties of the source object(s)
            
            nSrc = length(src);
            nC   = size(obj.c,2);
            for kSrc=1:nSrc
                src(kSrc).mask      = obj.pupilLogical;
                src(kSrc).amplitude = obj.pupil;
                switch nC
                    case 1
                        src(kSrc).phase     = utilities.toggleFrame(obj.p_p*obj.c*src(kSrc).waveNumber,3);
                    case nSrc
                        src(kSrc).phase     = utilities.toggleFrame(obj.p_p*obj.c(:,kSrc)*src(kSrc).waveNumber,3);
                    otherwise
                        src(kSrc).phase     = utilities.toggleFrame(obj.p_p*obj.c*src(kSrc).waveNumber,3);
                end
            end
        end
        
        function out = fourier(obj,f,o)
            %% FOURIER Zernike polynomials Fourier transform
            %
            % out = fourier(obj,f,o) computes the Fourier transform of the
            % Zernike polynomials at the spatial frequency vector [f,o]
            
            out = zeros(size(f,1),size(f,2),obj.nMode);
            for kMode = 1:obj.nMode
                out(:,:,kMode) = fun(obj.j,obj.n,obj.m);
            end
            function out1 = fun(zj,n,m)
                krkr = m~=0;
                g = (-1).*((n+m*krkr)/2).*1i.^(m.*krkr).*2.^(krkr/2);
                p = pi.*krkr.*((-1).^zj-1)/4;
                out1 = 2.*sqrt(obj.n+1).*utilities.sombrero(pi.*f.*obj.D).*...
                    g.*cos(m.*o+p);
            end
        end
        
        function out = fourierAzimSum(obj,f)
            %% FourierAzimSum azimuth summed Fourier transform
            %
            % out = fourierAzimSum(obj,f) computes the azimuth average of
            % the Fourier transform of the Zernike polynomials at the
            % spatial frequency f
            
            out = zeros(size(f,1),obj.nMode);
            for kMode = 1:obj.nMode
                out(:,kMode) = fun(obj.j,obj.n,obj.m);
            end
            function out1 = fun(zj,n,m)
                krkr = m~=0;
                g = (-1).*((n+m*krkr)/2).*1i.^(m.*krkr).*2.^(krkr/2);
                p = pi.*krkr.*((-1).^zj-1)/4;
                out1 = 2.*sqrt(obj.n+1).*utilities.sombrero(pi.*f.*obj.D).*g;
                if krkr
                    out1 = out1.*(sin(2.*pi.*m+p)-sin(p))./m;
                end
            end
        end
        
        function out = fourierPeakFreq(obj)
            %% FOURIERPEAKFREQ Frequency of Zernike maximum power
            %
            % out = fourierPeakFreq(obj) computes the frequency at which
            % the Zernike polynomials spectrum is maximal
            
            out = (obj.n+1)/(pi*obj.D);
        end
        
        %% get Zernike x derivative coefficients
        function out = get.xDerivativeMatrix(obj)
            i = []; j = []; s = [];
            nModes = length(obj.j);
            for kModes=2:nModes
                modesProjection = 1:kModes-1;
                rule  = abs( obj.m(kModes) - obj.m(modesProjection) ) == 1;
                mNnz = obj.m(kModes) ~= 0 & obj.m(modesProjection) ~= 0;
                rule1 = mNnz & ...
                    ~mod(obj.j(kModes) - obj.j(modesProjection),2);
                rule2 = (obj.m(kModes) == 0 & obj.m(modesProjection) ~= 0) & ...
                    ~rem(obj.j(modesProjection),2);
                rule3 = (obj.m(kModes) ~= 0 & obj.m(modesProjection) == 0) & ...
                    ~rem(obj.j(kModes),2);
                rule(~(rule1 | rule2 | rule3)) = 0;
                nnz = sum(rule);
                c   = sqrt( (obj.n(kModes)+1).*(obj.n(modesProjection(rule))+1) );
                c(~mNnz(rule)) = c(~mNnz(rule)).*sqrt(2);
                i   = [i modesProjection(rule)];
                j   = [j ones(1,nnz).*kModes];
                s   = [s c];
            end
            out = sparse(i,j,s,nModes,nModes);
        end
        
        %% get Zernike y derivative coefficients
        function out = get.yDerivativeMatrix(obj)
            i = []; j = []; s = [];
            nModes = length(obj.j);
            for kModes=2:nModes
                modesProjection = 1:kModes-1;
                rule  = abs( obj.m(kModes) - obj.m(modesProjection) ) == 1;
                mNnz = obj.m(kModes) ~= 0 & obj.m(modesProjection) ~= 0;
                rule1 = mNnz & ...
                    mod(obj.j(kModes) - obj.j(modesProjection),2);
                rule2 = (obj.m(kModes) == 0 & obj.m(modesProjection) ~= 0) & ...
                    rem(obj.j(modesProjection),2);
                rule3 = (obj.m(kModes) ~= 0 & obj.m(modesProjection) == 0) & ...
                    rem(obj.j(kModes),2);
                rule(~(rule1 | rule2 | rule3)) = 0;
                nnz = sum(rule);
                c   = sqrt( (obj.n(kModes)+1).*(obj.n(modesProjection(rule))+1) ); %#ok<*PROP>
                c(~mNnz(rule)) = c(~mNnz(rule)).*sqrt(2);
                signRule = obj.m(kModes) ~= 0 & ( ...
                    obj.m(modesProjection(rule)) == (obj.m(kModes)+1) & ...
                    rem(obj.j(kModes),2) );
                c(signRule) = -c(signRule);
                signRule = obj.m(kModes) ~= 0 & ( ...
                    obj.m(modesProjection(rule)) == (obj.m(kModes)-1) & ...
                    ~rem(obj.j(kModes),2) );
                c(signRule) = -c(signRule);
                i   = [i modesProjection(rule)];
                j   = [j ones(1,nnz).*kModes];
                s   = [s c];
            end
            out = sparse(i,j,s,nModes,nModes);
        end
        
        %% get Zernike x derivative
        function out = get.xDerivative(obj)
            out = obj.p_p*obj.xDerivativeMatrix;
            if ~obj.lex
                out = utilities.toggleFrame(out,3);
            end
        end
        
        %% get Zernike y derivative
        function out = get.yDerivative(obj)
            out = obj.p_p*obj.yDerivativeMatrix;
            if ~obj.lex
                out = utilities.toggleFrame(out,3);
            end
        end
        
        %% get noise covariance
        function out = get.noiseCovariance(obj)
            ggx = obj.xDerivativeMatrix;
            ggy = obj.yDerivativeMatrix;
            DtD = full(pi.*(ggx'*ggx + ggy'*ggy));
            out = pinv(DtD);
        end
        
        function out = circularCut(obj,delta,largeSmallRadiusRatio)
            xTrunc = obj.r.*cos(obj.o)/largeSmallRadiusRatio+delta(1);
            yTrunc = obj.r.*sin(obj.o)/largeSmallRadiusRatio+delta(2);
            rTrunc = hypot(xTrunc,yTrunc);
            oTrunc = atan2(yTrunc,xTrunc);
            rTrunc = rTrunc(:);
            oTrunc = oTrunc(:);
            nj = length(obj.j);
            out = zeros(length(rTrunc),nj);
            for kj = 1:nj
                out(:,kj) = zernike.fun(obj.j(kj),obj.n(kj),obj.m(kj),rTrunc,oTrunc).*(obj.r(:)<=1).*(rTrunc<=1);
            end
        end
        
        function P = smallFootprintExpansion(obj,delta,largeSmallRadiusRatio)
            %% SMALLFOOTPRINTEXPANSION Portion of Zernike expansion onto zernikes
            %
            % obj =
            % obj.smallFootprintExpansion(delta,largeSmallRadiusRatio)
            % expands a circular portion of Zernike polynomials centered on
            % delta(1) and delta(2) onto another Zernike basis fitted to
            % the circular portion. The ratio of the radius of both basis
            % is given by largeSmallRadiusRatio.
            
            nMode = length(obj.j);
            lastMode = max(obj.j);
            if all(delta==0) && largeSmallRadiusRatio==1
                P = diag(ones(1,nMode),nMode-lastMode);
                if lastMode>nMode
                    P(:,end+(nMode-lastMode)+1:end) = [];
                end
            else
                znmj = mat2cell([obj.j;obj.n;obj.m],3,ones(1,nMode));
                znmj = repmat(znmj,lastMode,1);
                zern = zernike(1:lastMode);
                znmi = mat2cell([zern.j;zern.n;zern.m],3,ones(1,lastMode));
                znmi = repmat(znmi,nMode,1);
                znmi = znmi';
                index = triu(true(lastMode,nMode),nMode-lastMode);
                pij = @(znmi,znmj) ...
                    quad2d(@(r,o) integrand(r,o,znmi(1),znmi(2),znmi(3),znmj(1),znmj(2),znmj(3)) ,0,1,0,2*pi);
                P = zeros(lastMode,nMode);
                %                 tic
                P(index) = cellfun(pij,znmj(index),znmi(index))./pi;
                %                 toc
                P(abs(P)<1e-6) = 0;
            end
            
            function out = integrand(r,o,zi,ni,mi,zj,nj,mj)
                xTrunc = r.*cos(o)/largeSmallRadiusRatio+delta(1);
                yTrunc = r.*sin(o)/largeSmallRadiusRatio+delta(2);
                rTrunc = hypot(xTrunc,yTrunc);
                oTrunc = atan2(yTrunc,xTrunc);
                out = r.*zernike.fun(zi,ni,mi,rTrunc,oTrunc) .* ...
                    zernike.fun(zj,nj,mj,r,o);
            end
            
        end
        
        function out = otf(obj, r)
            out = [];
        end
        function out = psf(obj,f)  
            out = [];
        end
        function out = fullWidthHalfMax(obj)
            out = [];
        end

        
    end
    
    % ---------------------------------------------------------------------
    methods (Static)
        
        function out = fun(j,n,m,r,o)
            %FUN Zernike polynomials mathematical expression
            %
            % out = zernike.fun(j,n,m,z) computes the Zernike polynomial
            % defined by the mode j, the radial order n, the azimuthal
            % frequency m at the locations given by the polar coordinate r
            % and o
            
            krkr = m~=0;
            out = sqrt(n+1).*R_fun(r,n,m).*2.^(0.5.*krkr).*...
                cos(m.*o+pi.*krkr.*((-1).^j-1)./4);
            
            %             out = sqrt(n+1)*R_fun(r,n,m);
            %             if m~=0
            %                 if rem(j,2)
            %                     out = out.*sqrt(2).*sin(m.*o);
            %                 else
            %                     out = out.*sqrt(2).*cos(m.*o);
            %                 end
            %             end
            
            function R = R_fun(r,n,m)
                R=zeros(size(r));
                for s=0:(n-m)/2
                    R = R + (-1).^s.*prod(1:(n-s)).*r.^(n-2.*s)./...
                        (prod(1:s).*prod(1:((n+m)/2-s)).*prod(1:((n-m)/2-s)));
                end
            end
            
        end
        
        function out = funOld(j,n,m,z)
            %FUN Zernike polynomials mathematical expression
            %
            % out = zernike.fun(j,n,m,z) computes the Zernike polynomial
            % defined by the mode j, the radial order n, the azimuthal
            % frequency m at the locations given by the polar coordinate r
            % and o
            
            r = abs(z);
            o = angle(z);
            
            out = sqrt(n+1)*R_fun(r,n,m);
            
            if m~=0
                if rem(j,2)
                    out = out.*sqrt(2).*sin(m.*o);
                else
                    out = out.*sqrt(2).*cos(m.*o);
                end
            end
            
            function R = R_fun(r,n,m)
                R=zeros(size(r));
                for s=0:(n-m)/2
                    R = R + (-1).^s.*prod(1:(n-s)).*r.^(n-2.*s)./...
                        (prod(1:s).*prod(1:((n+m)/2-s)).*prod(1:((n-m)/2-s)));
                end
            end
            
        end
        
        function out = nModeFromRadialOrder(n)
            % NMODEFROMRADIALORDER Number of Zernike polynomials
            %
            % out = zernike.nModeFromRadialOrder(n) returns the number of
            % Zernike polynomials (n+1)(n+2)/2 for a given radial order n
            out = (n+1).*(n+2)/2;
        end
        
        function P = demoProjection(delta,largeSmallRadiusRatio)
            firstMode = 2;
            lastMode = 15;
            nMode = lastMode - firstMode + 1;
            zernPlus  = zernike(firstMode:lastMode,32);
            zernMinus = zernike(1:lastMode,32);
            P = smallFootprintExpansion(zernPlus,delta,largeSmallRadiusRatio);
            size(P)
            zernPlusCut = circularCut(zernPlus,delta,largeSmallRadiusRatio);
            zernPlusCut = reshape(zernPlusCut,[32,32,nMode]);
            zernPlus.c = eye(nMode);
            zernPlus.lex = false;
            zernMinus.c = P;
            zernMinus.lex = false;
            figure
            h = imagesc([zernPlus.phase(:,:,1),...
                zernPlusCut(:,:,1),...
                zernMinus.phase(:,:,1)]);
            axis equal tight
            colorbar
            for kMode = 1:nMode
                set(h,'CData',[zernPlus.phase(:,:,kMode),...
                    zernPlusCut(:,:,kMode),...
                    zernMinus.phase(:,:,kMode)])
                title(num2str(kMode+firstMode-1))
                pause
            end
        end
        
    end
    
    % ---------------------------------------------------------------------
    methods (Access=private)
        
        function [nf,mf] = findNM(obj)
            % FINDNM Zernike radial and azimuthal order finder
            % [n,m] = findNM(mode) Computes the radial and azimuthal order of modes of
            % Zernike polynomials
            
            mode = obj.j;
            for counti=1:size(mode,1)
                for countj=1:size(mode,2)
                    
                    %ordre radial
                    n = 0;
                    %ordre azimuthal
                    m = 0;
                    
                    count = 0;
                    
                    while length(n)<mode(counti,countj)
                        
                        count = count + 1;
                        tmp = 0:count;
                        tmp = tmp(rem(count-tmp,2)==0);
                        
                        if all(tmp)
                            n = [n ones(1,2.*length(tmp)).*count];
                            tmp1 = tmp(ones(2,1),:);
                            m = [m reshape(tmp1,1,size(tmp1,1).*size(tmp1,2))];
                        else
                            n = [n ones(1,(2.*length(tmp))-1).*count];
                            tmp1 = tmp(ones(2,1),2:length(tmp));
                            m = [m 0 reshape(tmp1,1,size(tmp1,1).*size(tmp1,2))];
                        end
                        
                    end
                    
                    nf(counti,countj) = n(mode(counti,countj));
                    mf(counti,countj) = m(mode(counti,countj));
                    
                end
            end
            
            %             % From Malacara: Optical Shop testing
            %             mode = 1:max(obj.j);
            %             nf = ceil( ( -3 + sqrt( 9 + 8*( mode-1 ) ) )/2);
            %             mf = zeros(size(nf));
            %             k = 1;
            %             while k<=length(nf)
            %                 index = nf==nf(k);
            %                 mf(index) = sort( abs( nf(index) - 2*( obj.j(index) - (nf(index)+1).*nf(index)/2 -1 ) ) );
            %                 k = k + sum(index);
            %             end
            %             index = ismenber(mode,obj
        end
        
        
        function fun = polynomials(obj,unitNorm)
            % POLYNOMIALS Zernike polynomials
            % fun = polynomes(obj) Computes the Zernike polynomials for the
            % Zernike object  sampled on the polar coordinates arrays radius and
            % angle. The radius must be normalized to 1.
            
            % disp('Computing the polynomes ...')
            
            if isempty(obj.r) || isempty(obj.o)
                
                fun = [];
                
            else
                
                nv = obj.n;
                mv = obj.m;
                nf  = length(obj.j);
                fun = zeros(numel(obj.r),nf);
                r = obj.r(obj.pupilLogical);
                o = obj.o(obj.pupilLogical);
                
                % Null azimuthal order
                ind_m = find(mv==0);
                for cpt=ind_m
                    n = nv(cpt);
                    m = mv(cpt);
                    fun(obj.pupilLogical,cpt) = sqrt(n+1).*R_fun(r,n,m);
                end
                mod_mode = mod(obj.j,2);
                % Even polynomes
                ind_m = find(mod_mode==0 & mv~=0);
                for cpt=ind_m
                    n = nv(cpt);
                    m = mv(cpt);
                    fun(obj.pupilLogical,cpt) = sqrt(n+1).*R_fun(r,n,m).*sqrt(2).*cos(m.*o);
                end
                % Odd polynomes
                ind_m = find(mod_mode==1 & mv~=0);
                for cpt=ind_m
                    n = nv(cpt);
                    m = mv(cpt);
                    fun(obj.pupilLogical,cpt) = sqrt(n+1).*R_fun(r,n,m).*sqrt(2).*sin(m.*o);
                end
                
            end
            % disp(' ===>> DONE')
            
            if unitNorm
                fun = fun*diag(1./obj.nollNorm);
            end
            
            % Radial function
            function R = R_fun(r,n,m)
                R=zeros(size(r));
                for s=0:(n-m)/2
                    R = R + (-1).^s.*prod(1:(n-s)).*r.^(n-2.*s)./...
                        (prod(1:s).*prod(1:((n+m)/2-s)).*prod(1:((n-m)/2-s)));
                end
                
                
            end
            
        end
        
    end
    
end