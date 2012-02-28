classdef telescopeAbstract < handle
% Create a telescopeAbstract object
%
% tel = telescopeAbstract(D) creates a telescopeAbstract object from
% the diameter D.
%
% tel = telescopeAbstract(D,'parameter',value,...) creates a
% telescopeAbstract object from the diameter D and from optionnal
% parameter-value pair arguments. The parameters are obstructionRatio,
% fieldOfViewInArcsec, fieldOfViewInArcmin or resolution.
%
% This class should never be called directly. It is an abstract class.
% To define a telescope object, the telescope class should be used
% instead. Both the telescope and the zernike classes inherit from the
% telescopeAbstract class
%
% See also telescope and zernike

properties
    % diameter
    D;
    % central obstruction ratio
    obstructionRatio;
%         % conjugation altitude
%         conjugationHeight;
%         % focalisation distance
%         focalDistance;
    % elevation
    elevation;
    % field-of-view
    fieldOfView;
    % diameter resolution in pixel
    resolution;
    % phase listener
    phaseListener;
    % phase listener
    wavelengthListener;
    % wind shifted turbulence sampling time
    samplingTime;
    leap
end % properties


properties (Abstract)
    % tag
    tag;
end % properties(Abstract)


properties (Abstract , Dependent)% , SetAccess = private)
    % telescope pupil mask
    pupil;
end % properties (Abstract , Dependent)


properties (Dependent,SetAccess=private)
    % radius
    R;
    % telescope pupil mask
    pupilLogical;
%         % telescope area
%         area;
    % telescope area in pixels
    pixelArea;
end % properties (Dependent,SetAccess=private)


properties (Dependent)
    % optical aberrations seen by the telescope
    opticalAberration;
end % properties (Dependent)


properties (Access=protected)
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
    log;
    p_pupil;
    ppp;
    layerStep;
    phaseScreenWavelength;
end % properties (Access=protected)


methods


%% Constructor
function obj = telescopeAbstract(D,varargin)

  p = inputParser;
  p.addRequired('D', @isnumeric);
  p.addParamValue('obstructionRatio', 0, @isnumeric);
%             p.addParamValue('conjugationHeight', 0, @isnumeric);
%             p.addParamValue('focalDistance', Inf, @isnumeric);
  p.addParamValue('elevation', 0, @isnumeric);
  p.addParamValue('fieldOfViewInArcsec', [], @isnumeric);
  p.addParamValue('fieldOfViewInArcmin', [], @isnumeric);
  p.addParamValue('resolution', [], @isnumeric);
  p.addParamValue('samplingTime', [], @isnumeric);
  p.addParamValue('opticalAberration', [], @(x) true);
  p.parse(D, varargin{:});
  obj.D                = p.Results.D;
%             obj.conjugationHeight = p.Results.conjugationHeight;
%             obj.focalDistance = p.Results.focalDistance;
  obj.obstructionRatio = p.Results.obstructionRatio;
  obj.elevation = p.Results.elevation;
  if ~isempty(p.Results.fieldOfViewInArcsec)
      obj.fieldOfView      = p.Results.fieldOfViewInArcsec./cougarConstants.radian2arcsec;
  elseif ~isempty(p.Results.fieldOfViewInArcmin)
      obj.fieldOfView      = 60*p.Results.fieldOfViewInArcmin./cougarConstants.radian2arcsec;
  else
      obj.fieldOfView      = 0;
  end
  obj.resolution       = p.Results.resolution;
  obj.samplingTime = p.Results.samplingTime;
  obj.log = logBook.checkIn(obj);
  obj.opticalAberration = p.Results.opticalAberration;

end % telescopeAbstract()

%% Get the logical pupil mask
function pupilLogical = get.pupilLogical(obj)

  pupilLogical = logical(obj.pupil>0);

end % get.pupilLogical()

%% Get telescope radius
function out = get.R(obj)

  out = obj.D/2;

end % get.R()

%% Get telescope surface
function out = area(obj)
    out = pi*obj.R^2*(1-obj.obstructionRatio^2);
end

%% Get telescope surface in pixels
function out = get.pixelArea(obj)

  out = sum(obj.pupil(:));

end % get.pixelArea()

function out = diameterAt(obj,height)

  out = obj.D + 2.*height.*tan(obj.fieldOfView/2);

end % diameterAt(

function out = zernike(obj,modes,varargin)
%% ZERNIKE

  out = zernike(modes,obj.D,'resolution',obj.resolution, ...
                            'pupil',obj.pupil,varargin{:});

end % zernike()

function reset(obj)
%% RESET Reset the atmosphere phase screens
%
% reset(obj) resets the atmosphere random stream generator and
% re-computes the initial phase screens of the layers. These
% are the same than the phase screens computed during the
% atmosphere initialization process

  obj.atm.rngStream.State = obj.atm.initState;
  for kLayer=1:obj.atm.nLayer
    m_atm = slab(obj.atm,kLayer);
    fprintf('   Layer %d:\n',kLayer)
    fprintf('   -> Computing initial phase screen (D=%3.2fm,n=%dpx) ...', ...
            m_atm.layer.D,m_atm.layer.nPixel)
    obj.atm.layer(kLayer).phase = ...
    fourierPhaseScreen(m_atm,m_atm.layer.D,m_atm.layer.nPixel);
    fprintf('  Done \n')
  end

end % reset()

function draw(obj)
%% DRAW Reset the atmosphere phase screens
%
% draw(obj) reet the phase screens of the layers

  for kLayer=1:obj.atm.nLayer
    m_atm = slab(obj.atm,kLayer);
%                 fprintf('   Layer %d:\n',kLayer)
%                 fprintf('            -> Computing initial phase screen (D=%3.2fm,n=%dpx) ...',m_atm.layer.D,m_atm.layer.nPixel)
    obj.atm.layer(kLayer).phase = fourierPhaseScreen(m_atm,m_atm.layer.D,m_atm.layer.nPixel);
  %                 fprintf('  Done \n')
  end

end % draw()

function out = FT(obj,f)
%% FT Fourier transform
%
% out = FT(obj,f) computes the Fourier transform of the
% telescope pupil

  out   = ones(size(f)).*pi.*obj.D.^2.*(1-obj.obstructionRatio.^2)./4;
  index = f~=0;
  u = pi.*obj.D.*f(index);
  surface = pi.*obj.D.^2./4;
  out(index) = surface.*2.*besselj(1,u)./u;
  if obj.obstructionRatio>0
    u = pi.*obj.D.*obj.obstructionRatio.*f(index);
    surface = surface.*obj.obstructionRatio.^2;
    out(index) = out(index) - surface.*2.*besselj(1,u)./u;
  end
  out = out./(pi.*obj.D.^2.*(1-obj.obstructionRatio.^2)./4);

end % FT()

function out = entrappedEnergy(obj,eHalfSize,trap,psfOrOtf)
%% ENTRAPPEDENERGY Encircled of ensquared energy
%
% out = entrappedEnergy(obj,eHalfSize,trap) computes the
% entraped energy in a circle of radius eHalfSize if trap is
% set to 'circle' or in a square of half length eHalfSize if
% trap is set to 'square'

  if nargin<4
    psfOrOtf = 'psf';
  end
  switch lower(psfOrOtf)
    case 'otf'
      switch lower(trap)
        case 'circle'
          out = 2*pi*quadgk( ...
                @(r) r.*otf(obj,r).*...
                (2.*besselj(1,2*pi.*eHalfSize.*r)./ ...
                (2*pi.*eHalfSize.*r)),0,obj.D)*...
                pi*eHalfSize^2;
        case 'square'
          a = 2*eHalfSize;
          out = quad2d(...
                @(o,r) r.*otf(obj,r).*...
                (sin(pi.*r.*cos(o).*a)./(pi.*r.*cos(o).*a)).*...
                (sin(pi.*r.*sin(o).*a)./(pi.*r.*sin(o).*a)), ...
                0,2*pi,0,obj.D).*a.*a;
        otherwise
          error('cougar:telescope:entrapedEnergy',...
                'The trap is either a circle or a square!')
      end
   otherwise
     switch lower(trap)
       case 'circle'
         out = quadgk(@(x)x.*psf(obj,x),0,eHalfSize)*2*pi;
       case 'square'
%                             if isa(obj,'giantMagellanTelescope')
%                                 out = quad2d(@(x,y)psf(obj,hypot(x,y),atan2(y,x)),0,eHalfSize,0,eHalfSize)*4;
%                             else
          out = quad2d(@(x,y)psf(obj,hypot(x,y),atan2(y,x)),-eHalfSize, eHalfSize,-eHalfSize,eHalfSize);%*4;
%                             end
       otherwise
         error('cougar:telescope:entrapedEnergy',...
               'The trap is either a circle or a square!')
     end
  end

end % entrappedEnergy()

%% Set/Get for opticalAberration property
function set.opticalAberration(obj,val)

  obj.atm = val;
  if ~isempty(val) && isa(val,'atmosphere')
    obj.phaseListener = addlistener(obj.atm.layer(1),'phase','PostSet',...
                                    @(src,evnt) obj.imagesc );
    obj.phaseListener.Enabled = false;
    if ~isempty(obj.samplingTime)
      init(obj);
    end
    obj.leap = zeros(2,val.nLayer);
%                 obj.wavelengthListener = addlistener(obj.atm,'wavelength','PostSet',...
%                     @(src,evnt) obj.wavelengthScale );
  end

end % set.opticalAberration()

function out = get.opticalAberration(obj)

  out = obj.atm;

end % get.opticalAberration()

function varargout = update(obj)
    %% UPDATE Phase screens deplacement
    %
    % update(obj) moves the phase screens of each layer of one time
    % step in the direction of the corresponding wind vectors
    %
    % obj = update(obj) moves the phase screens and returns the
    % object

  if ~isempty(obj.atm) % uncorrelated phase screens

    if isinf(obj.samplingTime)

      for kLayer=1:obj.atm.nLayer
        obj.atm.layer(kLayer).phase = fourierPhaseScreen(slab(obj.atm,kLayer));
      end

    elseif ~(obj.atm.nLayer==1 && (obj.atm.layer.windSpeed==0 || ...
             isempty(obj.atm.layer.windSpeed) ) )
      for kLayer=1:obj.atm.nLayer

         % sampling in meter
        pixelLength = obj.atm.layer(kLayer).D./(obj.atm.layer(kLayer).nPixel-1);
        % 1 pixel increased phase sampling vector
        u0 = (-1:obj.atm.layer(kLayer).nPixel).*pixelLength;
        % phase sampling vector
        u = (0:obj.atm.layer(kLayer).nPixel-1).*pixelLength;

        % phase displacement in meter
        leap = [obj.windVx(kLayer) obj.windVy(kLayer)].* ...
                (obj.count(kLayer)+1).*obj.samplingTime;
        % phase displacement in pixel
        pixelLeap = leap/pixelLength;

        notDoneOnce = true;

        while any(pixelLeap>1) || notDoneOnce
          notDoneOnce = false;

          if obj.count(kLayer)==0
            % 1 pixel around phase increase
            Z = ...
            obj.atm.layer(kLayer).phase(obj.innerMask{kLayer}(2:end-1,2:end-1));
            X = ...
            obj.A{kLayer}*Z + ...
            obj.B{kLayer}*randn(obj.atm.rngStream,size(obj.B{kLayer},2),1);
            obj.mapShift{kLayer}(obj.outerMask{kLayer})  = X;
            obj.mapShift{kLayer}(~obj.outerMask{kLayer}) = ...
            obj.atm.layer(kLayer).phase(:);
          end

          % phase displacement (not more than 1 pixel)
          step   = min(abs(leap),pixelLength).*sign(leap);
%                             obj.layerStep(kLayer) = step;

          xShift = u - step(1);
          yShift = u - step(2);
%                             [xi,yi] = meshgrid(xShift,yShift);
          obj.atm.layer(kLayer).phase = ...
          spline2({u0,u0},obj.mapShift{kLayer},{yShift,xShift});

          leap = leap - step;
          pixelLeap = leap/pixelLength;

        end

        obj.count(kLayer) = rem(obj.count(kLayer)+1,obj.nShift(kLayer));

      end

    end

  end

  if nargout>0
    varargout{1} = obj;
  end

end % update()

function varargout = uplus(obj)
%% UPLUS + Update operator
%
% +obj updates the atmosphere phase screens
%
% obj = +obj returns the telescope object

  update(obj)
  if nargout>0
    varargout{1} = obj;
  end

end % uplus()

function obj = plus(obj,otherObj)
    %% + Add a component to the telescope
    %
    % obj = obj + otherObj adds an other object to the telescope
    % object

    obj.opticalAberration = otherObj;

end % plus()

function obj = minus(obj,otherObj)
%% - Remove a component from the telescope
%
% obj = obj - otherObj removes an other object to the telescope
% object

  if isa(obj.opticalAberration,class(otherObj))
  obj.opticalAberration = [];
  else
    warning('cougar:telescope:minus',...
    'The current and new objet must be from the same class (current: %s ~= new: %s)',...
    class(obj.opticalAberration),class(otherObj))
  end

end % minus()

function relay(obj,srcs)
%% RELAY Telescope to source relay
%
% relay(obj,srcs) writes the telescope amplitude and phase into
% the properties of the source object(s)

  if isempty(obj.resolution) % Check is resolution has been set
    if isscalar(srcs(1).amplitude) % if the src is not set either, do nothing
      return
    else % if the src is set, set the resolution according to src wave resolution
      obj.resolution = length(srcs(1).amplitude);
    end
  end

  nSrc = numel(srcs);
  for kSrc=1:nSrc % Browse the srcs array
    src = srcs(kSrc);
    % Set mask and pupil first
    src.mask      = obj.pupilLogical;
    if isempty(src.nPhoton) || (isempty(obj.samplingTime) || ...
       isinf(obj.samplingTime))
      src.amplitude = obj.pupil;
    else
      src.amplitude = ...
      obj.pupil.*sqrt(obj.samplingTime*src.nPhoton.*obj.area/sum(obj.pupil(:)));
    end
    out = 0;
    if ~isempty(obj.atm) % Set phase if an atmosphere is defined
      if obj.fieldOfView==0 && isNgs(src)
        out = out + sum(cat(3,obj.atm.layer.phase),3);
      else
        atm_m           = obj.atm;
        nLayer          = atm_m.nLayer;
        layers          = atm_m.layer;
        altitude_m      = [layers.altitude];
        sampler_m       = obj.sampler;
        phase_m         = { layers.phase };
        R_              = obj.R;
        layerSampling_m = obj.layerSampling;
        srcDirectionVector1 = src.directionVector(1);
        srcDirectionVector2 = src.directionVector(2);
        srcHeight = src.height;
        out = zeros(size(src.amplitude,1),size(src.amplitude,2),nLayer);
        parfor kLayer = 1:nLayer
        height = altitude_m(kLayer);
        [xs,ys] = meshgrid(layerSampling_m{kLayer});
        if height==0
          out(:,:,kLayer) = phase_m{kLayer};
        else
          layerR = R_*(1-height./srcHeight);
          u = sampler_m*layerR;
          xc = height.*srcDirectionVector1;
          yc = height.*srcDirectionVector2;
          [xi,yi] = meshgrid(u+xc,u+yc);
          out(:,:,kLayer) = linear(xs,ys,phase_m{kLayer},xi,yi);
        end
      end
      out = sum(out,3);
    end
      out = (obj.phaseScreenWavelength/src.wavelength)*out; % Scale the phase according to the src wavelength
    end
    src.phase = fresnelPropagation(src,obj) + out/sqrt( cos( obj.elevation ) );
    if isfinite(src.height);src.amplitude = 1./src.height;end
    src.timeStamp = src.timeStamp + obj.samplingTime;
  end

end % relay()

function relayReplacement(obj,srcs)

  if isempty(obj.resolution) % Check is resolution has been set
    if isscalar(srcs(1).amplitude) % if the src is not set either, do nothing
      return
    else % if the src is set, set the resolution according to src wave resolution
      obj.resolution = length(srcs(1).amplitude);
    end
  end

  nSrc   = numel(srcs);

  if isempty(obj.atm)
    for kSrc=1:nSrc % Browse the srcs array
      src = srcs(kSrc);
      src.mask      = obj.pupilLogical;
      if isempty(src.nPhoton)
        src.amplitude = obj.pupil;
      else
        src.amplitude = obj.pupil.*sqrt(obj.samplingTime*src.nPhoton.*obj.area/sum(obj.pupil(:)));
      end
    end
    out = 0;
  else

  atm_m           = obj.atm;
  nLayer          = atm_m.nLayer;
  layers          = atm_m.layer;
  altitude_m      = [layers.altitude];
  sampler_m       = obj.sampler;
  R_              = obj.R;

  for kSrc=1:nSrc % Browse the srcs array

    src = srcs(kSrc);
    src.mask      = obj.pupilLogical;
    if isempty(src.nPhoton)
      src.amplitude = obj.pupil;
    else
      src.amplitude = ...
      obj.pupil.* ...
      sqrt(obj.samplingTime*src.nPhoton.*obj.area/sum(obj.pupil(:)));
    end
    srcDirectionVector1 = src.directionVector(1);
    srcDirectionVector2 = src.directionVector(2);
    srcHeight = src.height;

    out = zeros(obj.resolution,obj.resolution,nLayer);
    for kLayer = 1:nLayer

      height = altitude_m(kLayer);
      % beam layer footprint diameter
      layerR = R_*(1-height./srcHeight);
      % and beam sampling
      u = sampler_m*layerR;
      % chief ray location
      xc = height.*srcDirectionVector1;
      yc = height.*srcDirectionVector2;
      % phase displacement in meter
      obj.leap(:,kLayer) = obj.leap(:,kLayer) + ...
          [obj.windVx(kLayer) ; obj.windVy(kLayer)].*obj.samplingTime;
      m_leap = obj.leap(:,kLayer);
      % beam mesh
      [xi,yi] = meshgrid(u+xc-m_leap(1),u+yc-m_leap(2));
      % pixel size
      p = layers(kLayer).D/(layers(kLayer).nPixel-1);
      % 2 pixels expanded layer phase screen radius
      bigR = p*(layers(kLayer).nPixel+1)/2;
      % and phase screen sampling
      [xs,ys] = meshgrid(linspace(-1,1,layers(kLayer).nPixel+2)*bigR);
      % check for out of range
      if any(abs(xi(:))>bigR | abs(yi(:))>bigR)
        % extract layer phase screen where to expand from
        u0 = 2:atm_m.layer(kLayer).nPixel+1;
        obj.atm.layer(kLayer).phase = obj.mapShift{kLayer}(u0-sign(m_leap(2)),u0-sign(m_leap(1)));
        % now 1 pixel expansion
        Z = obj.atm.layer(kLayer).phase(obj.innerMask{kLayer}(2:end-1,2:end-1));
        X = obj.A{kLayer}*Z + obj.B{kLayer}*randn(obj.atm.rngStream,size(obj.B{kLayer},2),1);
        obj.mapShift{kLayer}(obj.outerMask{kLayer})  = X;
        obj.mapShift{kLayer}(~obj.outerMask{kLayer}) = obj.atm.layer(kLayer).phase(:);
        % re-ajusting the phase displacement
        obj.leap(:,kLayer) = obj.leap(:,kLayer) - sign(m_leap).*p;
        m_leap = obj.leap(:,kLayer);
        [xi,yi] = meshgrid(u+xc-m_leap(1),u+yc-m_leap(2));
      end
      % phase interpolation
      out(:,:,kLayer) = linear(xs,ys,obj.mapShift{kLayer},xi,yi);

    end % layers
    out = sum(out,3);
    out = (obj.atm.wavelength/src.wavelength)*out; % Scale the phase according to the src wavelength

    end % srcs
  end % atm

  src.phase = fresnelPropagation(src,obj) + out;
  src.timeStamp = src.timeStamp + obj.samplingTime;

end % relayReplacement()

function varargout = imagesc(obj,varargin)
%% IMAGESC Phase screens display
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
    src = [];
    if nargin>1 && isa(varargin{1},'source')
      src = varargin{1};
      varargin(1) = [];
    end
    [n1,m1] = size(obj.atm.layer(1).phase);
    pupil = utilities.piston(n1,'type','logical');
    map = (obj.atm.layer(1).phase - mean(obj.atm.layer(1).phase(pupil))).*pupil;
    obj.imageHandle(1) = image([1,m1],[1,n1],map,...
        'CDataMApping','Scaled',varargin{:});
    hold on
    o = linspace(0,2*pi,101)';
    xP = obj.resolution*cos(o)/2;
    yP = obj.resolution*sin(o)/2;
    plot(xP+(n1+1)/2,yP+(n1+1)/2,'color',ones(1,3)*0.8)
    if ~isempty(src)
      kLayer = 1;
      for kSrc=1:numel(src)
        q = 1 - obj.atm.layer(kLayer).altitude/src(kSrc).height;
        xSrc = src(kSrc).directionVector(1).*...
               obj.atm.layer(kLayer).altitude.*...
               obj.atm.layer(kLayer).nPixel/...
               obj.atm.layer(kLayer).D;
        ySrc = src(kSrc).directionVector(2).*...
               obj.atm.layer(kLayer).altitude.*...
               obj.atm.layer(kLayer).nPixel/...
               obj.atm.layer(kLayer).D;
        plot(xSrc+xP*q+(n1+1)/2,ySrc+yP*q+(n1+1)/2,'color',ones(1,3)*0.8)
      end
    else
      plot(xP+(n1+1)/2,yP+(n1+1)/2,'k:')
    end
    text(m1/2,n1+0.5,...
         sprintf('%.1fkm: %.1f%%\n%.2fm - %dpx',...
         obj.atm.layer(1).altitude*1e-3,...
         obj.atm.layer(1).fractionnalR0*100,...
         obj.atm.layer(1).D,...
         obj.atm.layer(1).nPixel),...
         'HorizontalAlignment','Center',...
         'VerticalAlignment','Bottom')
    n = n1;
    offset = 0;
    for kLayer=2:obj.atm.nLayer
      [n,m] = size(obj.atm.layer(kLayer).phase);
      pupil = utilities.piston(n,'type','logical');
      offset = (n1-n)/2;
      map = (obj.atm.layer(kLayer).phase - ...
      mean(obj.atm.layer(kLayer).phase(pupil))).*pupil;
      obj.imageHandle(kLayer) = imagesc([1,m]+m1,[1+offset,n1-offset],map);
      if ~isempty(src)
        for kSrc=1:numel(src)
          q = 1 - obj.atm.layer(kLayer).altitude/src(kSrc).height;
          xSrc = src(kSrc).directionVector(1).*...
                 obj.atm.layer(kLayer).altitude.*...
                 obj.atm.layer(kLayer).nPixel/...
                 obj.atm.layer(kLayer).D;
          ySrc = src(kSrc).directionVector(2).*...
                 obj.atm.layer(kLayer).altitude.*...
                 obj.atm.layer(kLayer).nPixel/...
                 obj.atm.layer(kLayer).D;
          plot(xSrc+xP*q+m1+m/2,ySrc+yP*q+(n1+1)/2,'color',ones(1,3)*0.8)
        end
      else
        plot(xP+m1+m/2,yP+(n1+1)/2,'k:')
      end
      text(m1+m/2,(n1+1+m)/2,...
           sprintf('%.1fkm: %.1f%%\n%.2fm - %dpx',...
           obj.atm.layer(kLayer).altitude*1e-3,...
           obj.atm.layer(kLayer).fractionnalR0*100,...
           obj.atm.layer(kLayer).D,...
           obj.atm.layer(kLayer).nPixel),...
           'HorizontalAlignment','Center',...
           'VerticalAlignment','Bottom')
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

end % imagesc()

function wavelengthScale(obj,varargin)
%% WAVELENGTHSCALE

  if ~isempty(obj.mapShift)
    add(obj.log,obj,'Scaling wavefront!')
    nLayer = obj.atm.nLayer;
    for kLayer=1:nLayer
      obj.B{kLayer} = obj.atm.wavelengthScale.^2*obj.B{kLayer};
      obj.mapShift{kLayer} = obj.atm.wavelengthScale*obj.mapShift{kLayer};
    end
  end

end % wavelengthScale()


end % methods


methods (Abstract)
    display(obj)
    out = otf(obj, r)
    out = psf(obj,f)
    out = fullWidthHalfMax(obj)
end % methods (Abstract)


methods (Access=private)

function obj = init(obj)
%% INIT

  nInner = 2;
  obj.sampler = linspace(-1,1,obj.resolution);
  add(obj.log,obj,'Initializing phase screens making parameters:')
  obj.log.verbose = false;
  do = obj.D/(obj.resolution-1);
  for kLayer=1:obj.atm.nLayer
    if isempty(obj.atm.layer(kLayer).phase)
      D_m = obj.D + 2*obj.atm.layer(kLayer).altitude.*tan(0.5*obj.fieldOfView);
      nPixel = 1 + round(D_m./do);
      obj.atm.layer(kLayer).D = D_m;
      obj.atm.layer(kLayer).nPixel = nPixel;
      obj.layerSampling{kLayer}  = D_m*0.5*linspace(-1,1,nPixel);
      % ---------
      fprintf('   Layer %d:\n',kLayer)
      fprintf('   -> Computing initial phase screen (D=%3.2fm,n=%dpx) ...', ...
              D_m,nPixel)
      m_atm = slab(obj.atm,kLayer);
      obj.atm.layer(kLayer).phase = fourierPhaseScreen(m_atm,D_m,nPixel);
      fprintf('  Done \n')
      % ---------
      obj.outerMask{kLayer} = ...
          ~utilities.piston(nPixel,nPixel+2,...
          'shape','square','type','logical');
      obj.innerMask{kLayer} =  ...
          ~( obj.outerMask{kLayer} | ...
          utilities.piston(nPixel-2*nInner,nPixel+2,...
          'shape','square','type','logical') );
      fprintf('-> # of elements for the outer maks: %d and for the inner mask %d\n',...
      sum(obj.outerMask{kLayer}(:)),sum(obj.innerMask{kLayer}(:)));
      fprintf('-> Computing matrix A and B for layer %d: ',kLayer)
      [u,v] = meshgrid( (0:nPixel+1).*D_m/(nPixel-1) );
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
      pixelStep = [obj.windVx obj.windVy].*obj.samplingTime*(nPixel-1)/D_m;
      obj.nShift(kLayer) = max(floor(min(1./pixelStep)),1);
      u = (0:nPixel+1).*D_m./(nPixel-1);
      %                 [u,v] = meshgrid(u);
      obj.x{kLayer} = u;
      obj.y{kLayer} = u;%v;

      Z = obj.atm.layer(kLayer).phase(obj.innerMask{kLayer}(2:end-1,2:end-1));
      X = obj.A{kLayer}*Z + obj.B{kLayer}*randn(obj.atm.rngStream,size(obj.B{kLayer},2),1);
      obj.mapShift{kLayer}(obj.outerMask{kLayer})  = X;
      obj.mapShift{kLayer}(~obj.outerMask{kLayer}) = obj.atm.layer(kLayer).phase(:);

    end
  end
  obj.phaseScreenWavelength = obj.atm.wavelength;
  obj.log.verbose = true;

end % init()

end % methods( Access=Private )


methods (Static)


function out = symFT(symf)

  syms sD ri s
  u = pi*sD*symf;
  s = pi*sD^2/4;
  out = 2*s*besselj(1,u)/u;
  u = pi*sD*ri*symf;
  s = s*ri^2;
  out = out - 2*s*besselj(1,u)/u;
  out = out/(pi*sD^2*(1-ri^2)/4);

end % synFT()


end % methods( Static )


end classdef


function F = spline2(x,v,xi)
%2-D spline interpolation
% gridded spline interpolation via tensor products

  nv = size(v);
  d = length(x);
  values = v;
  sizeg = zeros(1,d);
  for i=d:-1:1
    ppp = spline(x{i},reshape(values,prod(nv(1:d-1)),nv(d)));
    values = ppval(ppp,xi{i}).';
    sizeg(i) = length(xi{i});
    nv = [sizeg(i), nv(1:d-1)];
  end
  F = reshape(values,sizeg);

end

function output = spline(x,y)

  output=[];

  sizey = size(y,1);
  n = length(x); yd = prod(sizey);

  % Generate the cubic spline interpolant in ppform

  dd = ones(yd,1); dx = diff(x); divdif = diff(y,[],2)./dx(dd,:);
  b=zeros(yd,n);
  b(:,2:n-1)=3*(dx(dd,2:n-1).*divdif(:,1:n-2)+dx(dd,1:n-2).*divdif(:,2:n-1));
  x31=x(3)-x(1);xn=x(n)-x(n-2);
  b(:,1)=((dx(1)+2*x31)*dx(2)*divdif(:,1)+dx(1)^2*divdif(:,2))/x31;
  b(:,n)= (dx(n-1)^2*divdif(:,n-2)+(2*xn+dx(n-1))*dx(n-2)*divdif(:,n-1))/xn;
  dxt = dx(:);
  c = spdiags([ [x31;dxt(1:n-2);0] ...
      [dxt(2);2*[dxt(2:n-1)+dxt(1:n-2)];dxt(n-2)] ...
      [0;dxt(2:n-1);xn] ],[-1 0 1],n,n);

  % sparse linear equation solution for the slopes
  mmdflag = spparms('autommd');
  spparms('autommd',0);
  s=b/c;
  spparms('autommd',mmdflag);
  % construct piecewise cubic Hermite interpolant
  % to values and computed slopes
  pp = pwch(x,y,s,dx,divdif); pp.dim = sizey;
  output = pp;

end % spline()


function F = linear(arg1,arg2,arg3,arg4,arg5)
%LINEAR 2-D bilinear data interpolation.
%   ZI = LINEAR(EXTRAPVAL,X,Y,Z,XI,YI) uses bilinear interpolation to
%   find ZI, the values of the underlying 2-D function in Z at the points
%   in matrices XI and YI.  Matrices X and Y specify the points at which
%   the data Z is given.  X and Y can also be vectors specifying the
%   abscissae for the matrix Z as for MESHGRID. In both cases, X
%   and Y must be equally spaced and monotonic.
%
%   Values of EXTRAPVAL are returned in ZI for values of XI and YI that are
%   outside of the range of X and Y.
%
%   If XI and YI are vectors, LINEAR returns vector ZI containing
%   the interpolated values at the corresponding points (XI,YI).
%
%   ZI = LINEAR(EXTRAPVAL,Z,XI,YI) assumes X = 1:N and Y = 1:M, where
%   [M,N] = SIZE(Z).
%
%   ZI = LINEAR(EXTRAPVAL,Z,NTIMES) returns the matrix Z expanded by
%   interleaving bilinear interpolates between every element, working
%   recursively for NTIMES. LINEAR(EXTRAPVAL,Z) is the same as
%   LINEAR(EXTRAPVAL,Z,1).
%
%   See also INTERP2, CUBIC.

  [nrows,ncols] = size(arg3);
  %     mx = numel(arg1); my = numel(arg2);
  s = 1 + (arg4-arg1(1))/(arg1(end)-arg1(1))*(ncols-1);
  t = 1 + (arg5-arg2(1))/(arg2(end)-arg2(1))*(nrows-1);

  % Matrix element indexing
  ndx = floor(t)+floor(s-1)*nrows;

  % Compute intepolation parameters, check for boundary value.
  if isempty(s), d = s; else d = find(s==ncols); end
  s(:) = (s - floor(s));
  if ~isempty(d), s(d) = s(d)+1; ndx(d) = ndx(d)-nrows; end

  % Compute intepolation parameters, check for boundary value.
  if isempty(t), d = t; else d = find(t==nrows); end
  t(:) = (t - floor(t));
  if ~isempty(d), t(d) = t(d)+1; ndx(d) = ndx(d)-1; end

  % Now interpolate.
  onemt = 1-t;
  F =  ( arg3(ndx).*(onemt) + arg3(ndx+1).*t ).*(1-s) + ...
       ( arg3(ndx+nrows).*(onemt) + arg3(ndx+(nrows+1)).*t ).*s;

end % linear()


function F = nearest(arg1,arg2,arg3,arg4,arg5)
%NEAREST 2-D Nearest neighbor interpolation.
%   ZI = NEAREST(EXTRAPVAL,X,Y,Z,XI,YI) uses nearest neighbor interpolation
%   to find ZI, the values of the underlying 2-D function in Z at the points
%   in matrices XI and YI.  Matrices X and Y specify the points at which
%   the data Z is given.  X and Y can also be vectors specifying the
%   abscissae for the matrix Z as for MESHGRID. In both cases, X
%   and Y must be equally spaced and monotonic.
%
%   Values of EXTRAPVAL are returned in ZI for values of XI and YI that are
%   outside of the range of X and Y.
%
%   If XI and YI are vectors, NEAREST returns vector ZI containing
%   the interpolated values at the corresponding points (XI,YI).
%
%   ZI = NEAREST(EXTRAPVAL,Z,XI,YI) assumes X = 1:N and Y = 1:M, where
%   [M,N] = SIZE(Z).
%
%   F = NEAREST(EXTRAPVAL,Z,NTIMES) returns the matrix Z expanded by
%   interleaving interpolates between every element.  NEAREST(EXTRAPVAL,Z)
%   is the same as NEAREST(EXTRAPVAL,Z,1).
%
%   See also INTERP2, LINEAR, CUBIC.

  [nrows,ncols] = size(arg3);
  mx = numel(arg1); my = numel(arg2);
  u = 1 + (arg4-arg1(1))/(arg1(mx)-arg1(1))*(ncols-1);
  v = 1 + (arg5-arg2(1))/(arg2(my)-arg2(1))*(nrows-1);

  % Check for out of range values of u and set to 1
  uout = (u<.5)|(u>=ncols+.5);
  anyuout = any(uout(:));
  if anyuout, u(uout) = 1; end

  % Check for out of range values of v and set to 1
  vout = (v<.5)|(v>=nrows+.5);
  anyvout = any(vout(:));
  if anyvout, v(vout) = 1; end

  % Interpolation parameters
  u = round(u); v = round(v);

  % Now interpolate
  ndx = v+(u-1)*nrows;
  F = arg3(ndx);

end % nearest()