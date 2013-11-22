classdef karhunenLoeve < hgsetget
    %% KARHUNENLOEVE Create a Karhunen-Loeve object
    %
    % obj = karhunenLoeve(hdfFile) creates a Karhunen-Loeve object from a
    % HDF file written by a CppAoApi program
    
    properties
        % azimuth vector
        angle;
        % azimuthal frequency compact vector
        azimOrder;
        % azimuthal frequency expanded vector
        m
        % lengh of azimuth vector
        na;
        % number of Karhunen-Loeve modes
        nf;
        % length of the radius vector
        nr;
        % number of radial function
        nv;
        % radial function matrix
        radialFun;
        % radius vector
        radius2;
        % variance of Karhunen-Loeve coefficients (compact vector)
        varCoef;
        % variance of Karhunen-Loeve coefficients (expanded vector)
        sigma2;
        % telescope diameter
        D;
        % telescope central obstruction ratio
        ri;
        % lexicographic ordering of frames (default: true)
        lex = true;
    end
    
    properties (Access=private)
        radialFunIndex;
    end
    
    methods
        
        %% Constructor
        function obj = karhunenLoeve(hdfFile)
            info = h5info(hdfFile);
            for k=1:length(info.Datasets)
                set(obj,info.Datasets(k).Name,...
                    h5read(hdfFile,['/',info.Datasets(k).Name]));
            end
            obj.D  = h5read(hdfFile,'/tel/D');
            obj.ri = h5read(hdfFile,'/tel/ri');
            obj.m = zeros(1,obj.nf);
            obj.sigma2 = zeros(1,obj.nf);
            obj.radialFunIndex = zeros(1,obj.nf);
            j = 0;
            for k=1:obj.nv
                j = j + 1;
                obj.radialFunIndex(j) = k;
                if obj.azimOrder(k)~=0
                    obj.m(j) = obj.azimOrder(k);
                    obj.sigma2(j) = obj.varCoef(k);
                    j = j + 1;
                    obj.radialFunIndex(j) = k;
                    obj.m(j) = obj.azimOrder(k);
                    obj.sigma2(j) = obj.varCoef(k);
                else
                    obj.sigma2(j) = obj.varCoef(k);
                end
            end
            obj.m(j) = obj.azimOrder(k);
            obj.nf = length(obj.m);
        end
        
        function out = pol2cart(obj,n_,j)
            %% POL2CART Cartesian grid interpolation
            %
            % modes = pol2cart(obj,n) interpolates the KL modes on a nxn
            % pixels cartesian grid
            %
            % modes = pol2cart(obj,n,j) interpolates the KL mode j on a nxn
            % pixels cartesian grid
            
            n = n_(1);    
            [x,y] = meshgrid( linspace(-1,1,n)*obj.D/2 );
            [o,r] = cart2pol(x,y);
            index = r>=obj.ri*obj.D/2 & r<=obj.D/2;
            r2 = r(index).^2;
            o = o(index);
            
            if nargin<3
                j = 1:obj.nf;
            end
            jv = j;
            nJ = length(jv);
            out = zeros(n^2,nJ);
            h = waitbar(0,'Computing the KL modes on a cartesian grid ...');
            for kJ=1:nJ
                j = jv(kJ);
                nRF = obj.radialFunIndex(j);
                radialFunInterp = interp1(obj.radius2,obj.radialFun(:,nRF),r2);
                m_m = obj.m(j);
                out(index,kJ) = radialFunInterp.*...
                    cos(m_m*o+pi.*(m_m~=0).*((-1).^j-1)./4);
                waitbar(kJ/nJ)
            end
            close(h)
            
            if length(n_)>1
                out = reshape(out,n_);
            end
            
        end
        
        function out = covarianceMatrix(obj,j)
            nRF = obj.radialFunIndex(j);
            var = obj.varCoef(nRF);
            nJ = length(j);
            out = spdiags(var,0,nJ,nJ);
        end
        
        
        function varargout = imagesc(obj,j)
            %% IMAGESC Display KL mode
            %
            % imagesc(obj,j) display KL mode j
            %
            % h = imagesc(...) return the image graphic handle
            
            n    = obj.radialFunIndex(j);
            m_m  = obj.m(j);
            mode = obj.radialFun(:,n)*...
                cos(obj.angle'*m_m + pi.*(m_m~=0).*((-1).^j-1)./4);
            [o,r] = meshgrid( obj.angle , sqrt(obj.radius2) );
            [x,y] = pol2cart(o,r);
            h = pcolor(x,y,mode);
            shading interp
            axis square
            colorbar
            if nargout>0
                varargout{1} = h;
            end
        end
        
    end
end