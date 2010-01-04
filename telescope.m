classdef telescope < telescopeCore
    % TELESCOPE Create a telescope object
    %
    % sys = telescope(D) creates a telescope object from the telescope diameter D
    %
    % sys = telescope(D,...) creates a telescope object from the
    % above parameter and from optionnal parameter-value pair arguments. The
    % optionnal parameters are those of the telescope clas.
    %
    % Example:
    % sys = telescope(10);
    % sys = telescope(30,'resolution',900)
    % tel = telescope(3.6,...
    %     'fieldOfViewInArcMin',2.5,...
    %     'resolution',nPx,...
    %     'samplingTime',1/100);
    %
    % See also telescopeCore
    
    properties
        % wind shifted turbulence sampling time
        samplingTime;
        % phase listener
        phaseListener;
    end
    
    properties (Dependent)
        % optical aberrations seen by the telescope
        opticalAberration;
    end
    
    properties (Access=private)
        atm;
        innerMask;
        outerMask;
        A;
        B;
        windVx;
        windVy;
        count;
        mapShift;
        nShift;
        x;
        y;
        imageHandle;
        layerSampling;
        sampler;
    end
    
    methods
        
        % Constructor
        function obj = telescope(D,varargin)
            p = inputParser;
            p.addRequired('D', @isnumeric);
            p.addParamValue('obstructionRatio', 0, @isnumeric);
            p.addParamValue('fieldOfViewInArcsec', [], @isnumeric);
            p.addParamValue('fieldOfViewInArcmin', [], @isnumeric);
            p.addParamValue('resolution', [], @isnumeric);
            p.addParamValue('samplingTime', [], @isnumeric);
            p.addParamValue('opticalAberration', [], @(x) isa(x,'atmosphere'));
            p.parse(D,varargin{:});
            obj = obj@telescopeCore(D,...
                'obstructionRatio',p.Results.obstructionRatio,...
                'fieldOfViewInArcsec',p.Results.fieldOfViewInArcsec,...
                'fieldOfViewInArcmin',p.Results.fieldOfViewInArcmin,...
                'resolution',p.Results.resolution);
            obj.samplingTime = p.Results.samplingTime;
            obj.opticalAberration = p.Results.opticalAberration;
        end
        
        % Destructor
        function delete(obj)
            if isa(obj.opticalAberration,'atmosphere')
                fprintf(' @(telescope) > Deleting atmosphere layer slabs!\n')
                for kLayer=1:obj.atm.nLayer
                    obj.atm.layer(kLayer).phase = [];
                end
            end
        end
        
        % Set/Get for opticalAberration property
        function set.opticalAberration(obj,val)
            obj.atm = val;
            if ~isempty(val) && isa(val,'atmosphere')
                obj.phaseListener = addlistener(obj.atm.layer(1),'phase','PostSet',...
                    @(src,evnt) obj.imagesc );
                obj.phaseListener.Enabled = false;
                init(obj);
            end
        end        
        function out = get.opticalAberration(obj)
            out = obj.atm;
        end
        
        function update(obj)
            % UPDATE Phase screens deplacement
            %
            % update(obj) moves the phase screens of each layer of one time
            % step in the direction of the corresponding wind vectors
            
%             disp(' (@telescope) > Layer translation')
            if ~isempty(obj.atm)
%                 disp('HERE')
                for kLayer=1:obj.atm.nLayer
                    
                    pixelLength = obj.atm.layer(kLayer).D./(obj.atm.layer(kLayer).nPixel-1); % sampling in meter
                    % 1 pixel increased phase sampling vector                    
                    u0 = (-1:obj.atm.layer(kLayer).nPixel).*pixelLength; 
                    % phase sampling vector
                    u = (0:obj.atm.layer(kLayer).nPixel-1).*pixelLength;

                    % phase displacement in meter
                    leap = [obj.windVx obj.windVy].*(obj.count(kLayer)+1).*obj.samplingTime;
                    % phase displacement im pixel
                    pixelLeap = leap/pixelLength;
                    
                    notDoneOnce = true;
                    while any(pixelLeap>1) || notDoneOnce
                        notDoneOnce = false;
                        
                        if obj.count(kLayer)==0
                            % 1 pixel around phase increase
                            Z = obj.atm.layer(kLayer).phase(obj.innerMask{kLayer}(2:end-1,2:end-1));
                            X = obj.A{kLayer}*Z + obj.B{kLayer}*randn(size(obj.B{kLayer},2),1);
                            obj.mapShift{kLayer}(obj.outerMask{kLayer})  = X;
                            obj.mapShift{kLayer}(~obj.outerMask{kLayer}) = obj.atm.layer(kLayer).phase(:);
                        end
                        
                        % phase displacement (not more than 1 pixel)
                        step   = min(leap,pixelLength);
                        xShift = u - step(1);
                        yShift = u - step(2);
                        obj.atm.layer(kLayer).phase ...
                               = spline2({u0,u0},obj.mapShift{kLayer},{yShift,xShift});
                        
                        leap = leap - step;
                        pixelLeap = leap/pixelLength;
                        
                    end
                        
                    obj.count(kLayer)       = rem(obj.count(kLayer)+1,obj.nShift(kLayer));
                    
                end
                
                
            end
            
        end
        
        function varargout = getPhaseScreen(obj,src)
            % GETPHASESCREEN Geometric phase propagation
            %
            % out = getPhaseScreen(obj,src) computes the geometric
            % propagation of the phase trough the turbulence layers in the
            % direction given by the source object src
            %
            % See also: source
            
            if nargin<2 || isempty(src) || isempty(obj.atm)
                out = zeros(size(obj.pupil));
            else
                if numel(src)>1
                    out = cell2mat( ...
                        cellfun( @obj.getPhaseScreen, num2cell(src), 'UniformOutput',false ) );
                else
                    if obj.fieldOfView==0 && isNgs(src)
                        out = sum(cat(3,obj.atm.layer.phase),3);
                    else
                        out = 0;
                        for kLayer = 1:obj.atm.nLayer
                            if obj.atm.layer(kLayer).altitude==0
                                out = out + obj.atm.layer(kLayer).phase;
                            else
                                layerR = obj.R*(1-obj.atm.layer(kLayer).altitude./src.height);
                                u = obj.sampler*layerR;
                                xc = obj.atm.layer(kLayer).altitude.*src.directionVector.x;
                                yc = obj.atm.layer(kLayer).altitude.*src.directionVector.y;
                                out = out + ...
                                    spline2({obj.layerSampling{kLayer},obj.layerSampling{kLayer}},...
                                    obj.atm.layer(kLayer).phase,{u+yc,u+xc});
                            end
                        end
                    end
                    out = out + fresnelPropagation(src,obj);
                    out(~obj.pupil) = 0;
                end
            end
            obj.pupilPhase = out;
            if nargout>0
                varargout{1} = out;
            end
        end
        
        function out = getWave(obj,src)
            if isempty(obj.atm) && numel(src)==0
                % if there is no optical aberration and no sources return
                % the telescope pupil
                out = ones(size(obj.pupil));
                out(~obj.pupil) = 0;
            elseif isempty(obj.atm)
                % if there is no optical aberration but at least 1 source:
                % fresnel wave propagation of the wave from the source to the
                % telescope pupil
                srcDims = size(src);
                out = bsxfun( @times , repmat( obj.pupil ,[1 srcDims(2)/srcDims(1)] ) ,...
                    exp( 1i.*cell2mat( ...
                    cellfun(@(x) fresnelPropagation(x,obj),num2cell(src),...
                    'UniformOutput',false) ) ) );
                if ~isempty([src.nPhoton])
                    photons = cell2mat( ...
                        cellfun( @(x) sqrt(x.nPhoton).*obj.area/sum(obj.pupil(:)),...
                        num2cell(src),'UniformOutput',false));
                    out = bsxfun( @times, out, photons);
                end
            elseif numel(src)>1
                % if there is optical aberration and more than 1 source
                out = cell2mat( ...
                    cellfun( @obj.getWave, num2cell(src), 'UniformOutput',false ) );
            else
                % if there is optical aberration and only 1 source
                out = exp(1i*getPhaseScreen(obj,src));
                out(~obj.pupil) = 0;
                if ~isempty(src.nPhoton)
                    out = out.*sqrt(src.nPhoton.*obj.area/sum(obj.pupil(:)));
                end
            end
        end
        
        function varargout = footprintProjection(obj,zernModeMax,src)
            nSource = length(src);
            P = cell(obj.atm.nLayer,nSource);
            for kSource = 1:nSource
                fprintf(' @(telescope) > Source #%2d - Layer #00',kSource)
                for kLayer = 1:obj.atm.nLayer
                    fprintf('\b\b%2d',kLayer)
                    obj.atm.layer(kLayer).zern = ...
                        zernike(1:zernModeMax,'resolution',obj.atm.layer(kLayer).nPixel);
                    conjD = obj.atm.layer(kLayer).D;
                    delta = obj.atm.layer(kLayer).altitude.*...
                        tan(src(kSource).zenith).*...
                        [cos(src(kSource).azimuth),sin(src(kSource).azimuth)];
                    delta = delta*2/conjD;
                    alpha = conjD./obj.D;
                    P{kLayer,kSource} = smallFootprintExpansion(obj.atm.layer(kLayer).zern,delta,alpha);
                    varargout{1} = P;
                end
                fprintf('\n')
            end
%             if nargout>1
%                 o = linspace(0,2*pi,101);
%                 varargout{2} = cos(o)./alpha + delta(1);
%                 varargout{3} = sin(o)./alpha + delta(2);
%             end
        end

        function varargout = imagesc(obj,varargin)
            % IMAGESC Phase screens display
            %
            % imagesc(obj) displays the phase screens of all the layers
            %
            % imagesc(obj,'PropertyName',PropertyValue) specifies property
            % name/property value pair for the image plot
            %
            % h = imagesc(obj,...) returns the image graphic handle
            
            if all(~isempty(obj.imageHandle)) && all(ishandle(obj.imageHandle))
                for kLayer=1:obj.atm.nLayer
                    n = size(obj.atm.layer(kLayer).phase,1);
                    pupil = utilities.piston(n,'type','logical');
                    map = (obj.atm.layer(kLayer).phase - mean(obj.atm.layer(kLayer).phase(pupil))).*pupil;
                    set(obj.imageHandle(kLayer),'Cdata',map);
                end
            else
                [n1,m1] = size(obj.atm.layer(1).phase);
                pupil = utilities.piston(n1,'type','logical');
                map = (obj.atm.layer(1).phase - mean(obj.atm.layer(1).phase(pupil))).*pupil;
                obj.imageHandle(1) = image([1,m1],[1,n1],map,...
                    'CDataMApping','Scaled',varargin{:});
                hold on
                o = linspace(0,2*pi,101);
                xP = obj.resolution*cos(o)/2;
                yP = obj.resolution*sin(o)/2;
                plot(xP+(n1+1)/2,yP+(n1+1)/2,'k:')
                n = n1;
                offset = 0;
                for kLayer=2:obj.atm.nLayer
                    [n,m] = size(obj.atm.layer(kLayer).phase);
                    pupil = utilities.piston(n,'type','logical');
                    offset = (n1-n)/2;
                    map = (obj.atm.layer(kLayer).phase - mean(obj.atm.layer(kLayer).phase(pupil))).*pupil;
                    obj.imageHandle(kLayer) = imagesc([1,m]+m1,[1+offset,n1-offset],map);
                    plot(xP+(m1+1)/2,yP+(n1+1)/2+offset/2,'k:')
                    m1 = m + m1;
                end
                hold off
                set(gca,'xlim',[1,m1],'ylim',[1+offset,n-offset],'visible','off')
                axis xy equal tight
                colorbar('location','southOutside')
            end
            if nargout>0
                varargout{1} = obj.imageHandle;
            end
        end
        
    end
    
    
    methods (Access=private)
        
        function obj = init(obj)
            nInner = 2;
            obj.sampler = linspace(-1,1,obj.resolution);
            for kLayer=1:obj.atm.nLayer
                if isempty(obj.atm.layer(kLayer).phase)
                    D = obj.D + 2*obj.atm.layer(kLayer).altitude.*tan(0.5*obj.fieldOfView);
                    obj.atm.layer(kLayer).D = D;
                    nPixel = ceil(1 + (obj.resolution-1)*D./obj.D);
                    obj.atm.layer(kLayer).nPixel = nPixel;
                    obj.layerSampling{kLayer}  = D*0.5*linspace(-1,1,nPixel);
                    % ---------
                    fprintf(' @(telescope)> Computing phase screen for layer %d (D=%3.2fm,n=%dpx) ...',kLayer,D,nPixel)
                    m_atm = slab(obj.atm,kLayer);
                    obj.atm.layer(kLayer).phase = fourierPhaseScreen(m_atm,D,nPixel);
                    fprintf('  Done \n')
                    % ---------
                    obj.outerMask{kLayer} = ...
                        ~utilities.piston(nPixel,nPixel+2,...
                        'shape','square','type','logical');
                    obj.innerMask{kLayer} =  ...
                        ~( obj.outerMask{kLayer} | ...
                        utilities.piston(nPixel-2*nInner,nPixel+2,...
                        'shape','square','type','logical') );
                    fprintf(' (@telescope)> # of elements for the outer maks: %d and for the inner mask %d\n',...
                        sum(obj.outerMask{kLayer}(:)),sum(obj.innerMask{kLayer}(:)));
                    fprintf(' @(telescope)> Computing matrix A and B for layer %d: ',kLayer)
                    [u,v] = meshgrid( (0:nPixel+1).*D/(nPixel-1) );
                    % ---------
                    innerZ = complex(u(obj.innerMask{kLayer}),v(obj.innerMask{kLayer}));
                    fprintf('ZZt ...')
                    ZZt = phaseStats.covarianceMatrix(innerZ,m_atm);
                    % ---------
                    outerZ = complex(u(obj.outerMask{kLayer}),v(obj.outerMask{kLayer}));
                    fprintf('\b\b\b, ZXt ...')
                    ZXt = phaseStats.covarianceMatrix(innerZ,outerZ,m_atm);
                    clear innerZ
                    % ---------
                    obj.A{kLayer}   = ZXt'/ZZt;
                    % ---------
                    clear ZZt
                    fprintf('\b\b\b, XXt ...')
                    XXt = phaseStats.covarianceMatrix(outerZ,m_atm);
                    clear outerZ
                    % ---------
                    BBt = XXt - obj.A{kLayer}*ZXt;
                    clear XXt ZXt
                    obj.B{kLayer} = chol(BBt,'lower');
                    fprintf('  Done \n')
                    % ---------
                    obj.windVx(kLayer) = m_atm.layer.windSpeed.*cos(m_atm.layer.windDirection);
                    obj.windVy(kLayer) = m_atm.layer.windSpeed.*sin(m_atm.layer.windDirection);
                    obj.count(kLayer) = 0;
                    obj.mapShift{kLayer} = zeros(nPixel+2);
                    pixelStep = [obj.windVx obj.windVy].*obj.samplingTime*(nPixel-1)/D;
                    obj.nShift(kLayer) = max(floor(min(1./pixelStep)),1);
                    u = (0:nPixel+1).*D./(nPixel-1);
                    %                 [u,v] = meshgrid(u);
                    obj.x{kLayer} = u;
                    obj.y{kLayer} = u;%v;
                end
            end
        end
        
    end
    
    methods (Static)
        function sys = demo(action)
            sys = telescope(2.2e-6,0.8,30,25,...
                'altitude',[0,10,15].*1e3,...
                'fractionnalR0',[0.7,0.2,0.1],...
                'windSpeed',[10,5,15],...
                'windDirection',[0,pi/4,pi/2],...
                'fieldOfViewInArcMin',2,...
                'resolution',60*15,...
                'samplingTime',1/500);
            %             sys = telescope(2.2e-6,0.8,10,25,...
            %                 'altitude',10.*1e3,...
            %                 'fractionnalR0',1,...
            %                 'windSpeed',10,...
            %                 'windDirection',0,...
            %                 'fieldOfViewInArcMin',1,...
            %                 'resolution',60,...
            %                 'samplingTime',1e-3);
            if nargin>0
                update(sys);
                sys.phaseListener.Enabled = true;
                imagesc(sys)
                while true
                    update(sys);
                    drawnow
                end
            end
        end
        
        function sys = demoSingleLayer
            sys = telescope(2.2e-6,0.8,10,25,...
                'altitude',10e3,...
                'fractionnalR0',1,...
                'windSpeed',10,...
                'windDirection',0,...
                'fieldOfViewInArcMin',2,...
                'resolution',256,...
                'samplingTime',1/500);
        end
    end
    
end

function F = spline2(x,v,xi)
%2-D spline interpolation

% % Determine abscissa vectors
% varargin{1} = varargin{1}(1,:);
% varargin{2} = varargin{2}(:,1).';
%
% %
% % Check for plaid data.
% %
% xi = varargin{4}; yi = varargin{5};
% xxi = xi(1,:); yyi = yi(:,1);
%
% %     F = splncore(varargin(2:-1:1),varargin{3},{yyi(:).' xxi},'gridded');
%     x = varargin(2:-1:1);
%     v = varargin{3};
%     xi = {yyi(:).' xxi};
% gridded spline interpolation via tensor products
nv = size(v);
d = length(x);
values = v;
sizeg = zeros(1,d);
for i=d:-1:1
    values = spline(x{i},reshape(values,prod(nv(1:d-1)),nv(d)),xi{i}).';
    sizeg(i) = length(xi{i});
    nv = [sizeg(i), nv(1:d-1)];
end
F = reshape(values,sizeg);

    function output = spline(x,y,xx)
        % disp('Part 1')
        % tStart = tic;
        output=[];
        
        sizey = size(y,1);
        n = length(x); yd = prod(sizey);
        
        % Generate the cubic spline interpolant in ppform
        
        dd = ones(yd,1); dx = diff(x); divdif = diff(y,[],2)./dx(dd,:);
        b=zeros(yd,n);
        b(:,2:n-1)=3*(dx(dd,2:n-1).*divdif(:,1:n-2)+dx(dd,1:n-2).*divdif(:,2:n-1));
        x31=x(3)-x(1);xn=x(n)-x(n-2);
        b(:,1)=((dx(1)+2*x31)*dx(2)*divdif(:,1)+dx(1)^2*divdif(:,2))/x31;
        b(:,n)=...
            (dx(n-1)^2*divdif(:,n-2)+(2*xn+dx(n-1))*dx(n-2)*divdif(:,n-1))/xn;
        dxt = dx(:);
        c = spdiags([ [x31;dxt(1:n-2);0] ...
            [dxt(2);2*[dxt(2:n-1)+dxt(1:n-2)];dxt(n-2)] ...
            [0;dxt(2:n-1);xn] ],[-1 0 1],n,n);
        
        % sparse linear equation solution for the slopes
        mmdflag = spparms('autommd');
        spparms('autommd',0);
        s=b/c;
        spparms('autommd',mmdflag);
        % toc(tStart)
        % construct piecewise cubic Hermite interpolant
        % to values and computed slopes
        %    disp('Part pwch')
        %    tStart = tic;
        pp = pwch(x,y,s,dx,divdif); pp.dim = sizey;
        % toc(tStart)
        
        % end
        
        %      disp('Part ppval')
        %  tStart = tic;
        output = ppval(pp,xx);
        % toc(tStart)
        
        
        
    end

end