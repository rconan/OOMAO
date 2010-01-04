classdef karhunenLoeve
    
    properties
        %         % radial order
        %         n;
        % azimuthal frequency
        m;
        % normalized radius
        r;
        % angle
        o;
        % modes radial functions
        pr;
        % coefficients
        c;
        % variance of the radial functions
        radialVariance;
        % lexicographic ordering of frames (default: true)
        lex = true;
    end
    
    properties (Dependent)
        % number of modes
        nMode;
        % radial function resolution
        nr;
        % azimuthal function resolution
        no;
        % azimuthal frequency resolution
        nm;
        %         % mode number
        %         j;
        % x coordinates of polar grid
        x
        % y coordinates of polar grid
        y
        % modes (in polar coordinates)
        p;
        % modes variance;
        variance;
    end
    
    methods
        
        % Constructor
        function obj = karhunenLoeve(m,atm,tel)
            error(nargchk(2,3,nargin));
            if nargin<3
                diameter  = atm.D;
                nPx = atm.resolution;
            else
                diameter  = tel.D;
                nPx = tel.resolution;
            end
            
            %             m = round(3*lastJ/4);
            nThetaPx = 8*m;
            theta        = zeros(1,1,nThetaPx);
            theta(1,1,:) = linspace(0,2*pi,8*m);
            r = sqrt(linspace(0,1,nPx)*diameter^2/4);
            % Kernel
            disp(' (@karhunenLoeve)> Kernel')
            r2_p_rp2 = bsxfun(@plus,r.^2,r'.^2);
            r_t_rp   = bsxfun(@times,r,r');
            rho = sqrt( bsxfun(@minus,r2_p_rp2,bsxfun(@times,r_t_rp,2*cos(theta))) );
            Gamma = phaseStats.covariance(rho,atm);
            clear r2_p_rp2 r_t_rp rho
            kernel = 0.5*real(fft(Gamma,[],3))./sqrt(nPx);
            % kernel = bsxfun(@times,r',kernel);
            % Eigen expansion
            eigenModes  = zeros(nPx,nPx,m);
            eigenValues = zeros(nPx,m);
            disp(' (@karhunenLoeve)> Eigen values derivation')
            for k=1:m
                [V,D] = eig(kernel(:,:,k));
                eigenModes(:,:,k) = V;
                eigenValues(:,k)  = diag(D);
            end
            % Modes builder
            [sortedEigenValues,indexEigenValues] = sort(eigenValues(:),'descend');
            eigenModes = reshape(eigenModes,[nPx,nPx*m]);
            area = pi*diameter^2/4;
            radialModes = eigenModes(:,indexEigenValues(1:m))/sqrt(area);
            %             radialOrder = repmat((nPx-1:-1:0)',1,nThetaPx);
            %             radialOrder = radialOrder(indexEigenValues(1:m));
            azimuthalFrequency = repmat(0:nThetaPx-1,nPx,1);
            azimuthalFrequency = azimuthalFrequency(indexEigenValues(1:m));
            theta = squeeze(theta);
            
            %             obj.n = radialOrder;
            obj.m = azimuthalFrequency;
            obj.r = r; %#ok<*PROP>
            obj.o = theta;
            obj.pr = radialModes;
            obj.radialVariance = eigenValues(indexEigenValues(1:m));
            %             obj.c = ones(length(obj.m),1);
        end
        
        function out = get.nMode(obj)
            out = 2*sum(obj.m~=0) + sum(obj.m==0);
        end
        
        function out = get.nr(obj)
            out = numel(obj.r);
        end
        
        function out = get.no(obj)
            out = numel(obj.o);
        end
        
        function out = get.nm(obj)
            out = numel(obj.m);
        end
        
        function out = get.x(obj)
            out = obj.r'*cos(obj.o');
        end
        
        function out = get.y(obj)
            out = obj.r'*sin(obj.o');
        end
        
        function basis = get.p(obj)
            azimuthalModes = bsxfun(@times,obj.o,obj.m')';
            basis = zeros(obj.nr,obj.no,obj.nMode);
            kModes = 0;
            for k=1:obj.nm
                kModes = kModes+1;
                if obj.m(k)==0;
                    basis(:,:,kModes) = bsxfun(@times,obj.pr(:,k),ones(1,obj.no));
                else
                    basis(:,:,kModes)   = bsxfun(@times,obj.pr(:,k),sqrt(2)*cos(azimuthalModes(k,:)));
                    kModes = kModes+1;
                    basis(:,:,kModes)   = bsxfun(@times,obj.pr(:,k),sqrt(2)*sin(azimuthalModes(k,:)));
                end
            end
        end
        
        function out = get.variance(obj)
            out = obj.radialVariance';
            out(obj.m==0) = NaN;
            out = [obj.radialVariance';out];
            out = out(:);
            out(isnan(out)) = [];
        end
        
        function basisCart = pCart(obj)
            
            nPx = obj.nr;
            [rc,oc] = utilities.cartAndPol(nPx,obj.r(end),'output','polar');
            index = rc<=obj.r(end);
            rc = rc(index);
            oc = oc(index);
            
            basisCart = zeros(nPx*nPx,obj.nMode);
            
            kModes = 0;
            for k=1:obj.nm
                kModes = kModes+1;
                buffer = spline(obj.r,obj.pr(:,k),rc);
                % radialModes norm must be 1 both in cartesian and polar coordinates
                buffer = buffer./norm(buffer);
                if obj.m(k)~=0;
                    radialModesCart = buffer;
                    buffer = radialModesCart.*sqrt(2).*cos(obj.m(k).*oc);
                    basisCart(index,kModes) = buffer;
                    kModes = kModes+1;
                    buffer = radialModesCart.*sqrt(2).*sin(obj.m(k).*oc);
                end
                basisCart(index,kModes) = buffer;
            end
        end
    end
end