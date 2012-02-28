classdef atmosphere < hgsetget
% Create an atmosphere object
%
% atm = atmosphere(wavelength,r0) creates an atmosphere object from the
% wavelength and the Fried parameter r0
%
% atm = atmosphere(wavelength,r0,L0) creates an atmosphere object from the
% wavelength, the Fried parameter r0 and the outer scale L0
%
% atm =
% atmosphere(wavelength,r0,'altitude',altitude,'fractionnalR0',fractionnalR
% 0) creates an atmosphere object from the wavelength, the Fried parameter
% r0, and from the altitudes and the fractionnalR0 of the turbulence layers
%
% atm =
% atmosphere(wavelength,r0,'altitude',altitude,'fractionnalR0',fractionnalR
% 0,'windSpeed',windSpeed,'windDirection',windDirection) creates an
% atmosphere object from the wavelength, the Fried parameter r0, and from
% the altitudes, the fractionnalR0, the wind speeds and the wind directions
% of the turbulence layers
%
% atm =
% atmosphere(wavelength,r0,L0,'altitude',altitude,'fractionnalR0',fractionn
% alR0) creates an atmosphere object from the wavelength, the Fried
% parameter r0, the outer scale L0 and from the altitudes and the
% fractionnalR0 of the turbulence layers
%
% atm =
% atmosphere(wavelength,r0,L0,'altitude',altitude,'fractionnalR0',fractionn
% alR0,'windSpeed',windSpeed,'windDirection',windDirection) creates an
% atmosphere object from the wavelength, the Fried parameter r0, the outer
% scale L0 and from the altitudes, the fractionnalR0, the wind speeds and
% the wind directions of the turbulence layers
%
% Example:
%     atm = atmosphere(photometry.V,0.15,30,...
%     'altitude',4e3,...
%     'fractionnalR0',1,...
%     'windSpeed',15,...
%     'windDirection',0);
%
% A full volume of turbulence above a telescope is made combining both
% a telescope and an atmosphere class
%
% The geometric propagation of a wavefront throught the atmosphere is
% obtained with a source object
%
% See also telescope and source


properties
    % Fried parameter
    r0;
    % turbulence outer scale
    L0;
    % number of turbulence layers
    nLayer;
    % turbulence layer object array
    layer;
    % atmosphere tag
    tag = 'ATMOSPHERE';
    % coherence function decay for theta0 and tau0 computation
    % . Roddier (default): coherence function 1/e decay
    % . Fried: structure function equal 1rd^2 (1/\sqrt(e) decay)
    coherenceFunctionDecay = exp(-1);
    % random number generator stream
    rngStream;
    % state at initialization of the atmosphere random number generator stream
    initState;
    wavelengthScale;
end % properties


properties (Dependent, SetObservable=true)
    % wavelength of the r0
    wavelength;
end % properties (Dependent, SetObservable=true)


properties (Dependent, SetAccess=private)
    % seeing
    seeingInArcsec;
    % theta0
    theta0InArcsec;
    % tau0
    tau0InMs;
    % mean height
    meanHeight;
    % mean velocity
    meanWind;
    % Greenwood frequency
    greenwoodFrequency;
end % properties (Dependent, SetAccess=private)


properties (Access=private)
    p_wavelength;
    p_log;
end % properties (Access=private)


methods


%% Constructor
function obj = atmosphere(lambda,r0,varargin)

  p = inputParser;
  p.addRequired('wavelength', @(x) isnumeric(x) || isa(x,'photometry') );
  p.addRequired('r0', @isnumeric);
  p.addOptional('L0', Inf, @isnumeric);
  p.addParamValue('altitude', 0, @isnumeric);
  p.addParamValue('fractionnalR0', 1, @isnumeric);
  p.addParamValue('windSpeed', [], @isnumeric);
  p.addParamValue('windDirection', [], @isnumeric);
  p.addParamValue('logging', true, @islogical);
  p.addParamValue('randStream', [], @(x) isempty(x) || isa(x,'RandStream') );
  p.parse(lambda,r0, varargin{:});
  if isa(p.Results.wavelength,'photometry')
    obj.p_wavelength = p.Results.wavelength.wavelength;
  else
    obj.p_wavelength = p.Results.wavelength;
  end
  obj.r0 = p.Results.r0;
  obj.L0 = p.Results.L0;
  obj.nLayer = length(p.Results.altitude);
  if any(isempty(p.Results.windSpeed))
    obj.layer = turbulenceLayer(p.Results.altitude, p.Results.fractionnalR0);
  else
    obj.layer = turbulenceLayer(...
                                p.Results.altitude,...
                                p.Results.fractionnalR0,...
                                p.Results.windSpeed,...
                                p.Results.windDirection);
  end
  if p.Results.logging
    obj.p_log = logBook.checkIn(obj);
    display(obj);
  end
  if isempty(p.Results.randStream)
    obj.rngStream = RandStream('mt19937ar');
  else
    obj.rngStream = p.Results.randStream;
  end
  obj.initState = obj.rngStream.State;

end % atmosphere()

% Destructor
function delete(obj)

  if ~isempty(obj.p_log)
      checkOut(obj.p_log,obj)
  end

end % delete()

function newObj = slab(obj,layerIndex)
%% SLAB Create a single turbulence layer atmosphere object
%
% singledAtm = slab(atm,k) creates an atmosphere object from
% the old atm object and the k-th turbulent layer

  newObj = atmosphere(...
                      obj.wavelength,...
                      obj.r0,...
                      obj.L0,...
                      'randStream',obj.rngStream,...
                      'logging',false);
  newObj.layer = obj.layer(layerIndex);

end % slab()

function display(obj)
%% DISPLAY Display object information
%
% disp(obj) prints information about the atmosphere object

  fprintf('___ %s ___\n',obj.tag)
  if isinf(obj.L0)
    fprintf(' Kolmogorov-Tatarski atmospheric turbulence:\n')
    fprintf('  . wavelength  = %5.2fmicron,\n  . r0          = %5.2fcm,\n  . seeing      = %5.2farcsec,\n  ',...
            obj.wavelength*1e6,obj.r0*1e2,obj.seeingInArcsec)
  else
    fprintf(' Von Karman atmospheric turbulence\n')
    fprintf('  . wavelength  = %5.2fmicron,\n  . r0          = %5.2fcm,\n  . L0          = %5.2fm,\n  . seeing      = %5.2farcsec,\n  ',...
            obj.wavelength*1e6,obj.r0*1e2,obj.L0,obj.seeingInArcsec)
  end
  if ~isinf(obj.theta0InArcsec)
    fprintf('. theta0(%2.0f%%) = %5.2farcsec,\n  ',...
    100*obj.coherenceFunctionDecay,obj.theta0InArcsec)
  end
  if ~isempty(obj.tau0InMs)
    fprintf('. tau0(%2.0f%%)   = %5.2fmillisec',...
            100*obj.coherenceFunctionDecay,obj.tau0InMs)
  else
    %fprintf('\b\b')
  end
  fprintf('\n')
  fprintf('----------------------------------------------------\n')
  fprintf('  Layer   Altitude[m]   fr0    wind([m/s] [deg])   D[m]    res[px]\n')
  for kLayer=1:obj.nLayer
    fprintf('  %2d      %8.2f      %4.2f    (%5.2f %6.2f)     %5.2f    %3d\n',...
            kLayer,...
            obj.layer(kLayer).altitude,...
            obj.layer(kLayer).fractionnalR0,...
            obj.layer(kLayer).windSpeed,...
            obj.layer(kLayer).windDirection*180/pi,...
            obj.layer(kLayer).D,...
            obj.layer(kLayer).nPixel)
  end
  fprintf('----------------------------------------------------\n')

end % display()

%% Set/Get wavelength property
function val = get.wavelength(obj)

  val = obj.p_wavelength;

end % get.wavelength(

function set.wavelength(obj,val)

  if isa(val,'photometry')
      val = val.wavelength;
  end
  obj.wavelengthScale = obj.wavelength/val;
  obj.r0 = obj.r0.*(val/obj.wavelength)^1.2;
  obj.p_wavelength = val;

end % set.wavelength()

%% Get seeingInArcsec property
function out = get.seeingInArcsec(obj)

  out = cougarConstants.radian2arcsec.*...
        0.98.*obj.wavelength./obj.r0;

end % get.seeingInArcsec()

%% Set coherenceFunctionDecay property
function set.coherenceFunctionDecay(obj,val)

  switch lower(val)
    case 'roddier'
      obj.coherenceFunctionDecay = exp(-1);
    case 'fried'
      obj.coherenceFunctionDecay = exp(-0.5);
    otherwise
      if isnumeric(val)
        obj.coherenceFunctionDecay = val;
      else
        error('oomao:atmosphere:coherenceFunctionDecay',...
              'The valid definitions for theta0 and tau0 are either', ...
              'Roddier, Fried or a numeric value between 0 and 1 excluded.')
      end
  end

end % set.coherenceFunctionDecay()

%% Get theta0InArcsec property
function out = get.theta0InArcsec(obj)

  z = [obj.layer.altitude];
  if (numel(z)==1 && z==0) || all([obj.layer.altitude]==0)
    out  = Inf;
  else
    if isinf(obj.L0)
      cst = -log(obj.coherenceFunctionDecay)*(24*gamma(6/5)/5)^(-5/6).* ...
            obj.r0.^(5./3);
      out = ( cst./sum( [obj.layer.fractionnalR0].*z.^(5./3) ) ).^(3/5);
    else
      fun = @(x) phaseStats.angularStructureFunction(abs(x),obj) + ...
            2.*log(obj.coherenceFunctionDecay);
      out = abs(fzero(fun,0));
    end
  end
  out = out*cougarConstants.radian2arcsec;

end % get.theta0InArcsec()

%% Get tau0InMs property
function out = get.tau0InMs(obj)

  v = [obj.layer.windSpeed];
  if isempty(v)
    out = [];
  elseif numel(v)==1 && v==0
    out  = Inf;
  else
    if isinf(obj.L0)
      cst = -log(obj.coherenceFunctionDecay)*(24*gamma(6/5)/5)^(-5/6).* ...
            obj.r0.^(5./3);
      out= (cst./sum( [obj.layer.fractionnalR0].*v.^(5./3) ) ).^(3/5);
    else
      fun = @(x) phaseStats.temporalStructureFunction(abs(x),obj) + ...
            2.*log(obj.coherenceFunctionDecay);
      out = abs(fzero(fun,0));
    end
  end
  out = out*1e3;

end % get.tau0InMs()

%% Get mean height
function out = get.meanHeight(obj)

  z = [obj.layer.altitude];
  out = sum( [obj.layer.fractionnalR0].*z.^(5./3) )^(3/5);

end % get.meanHeight()

%% Get mean wind velocity
function out = get.meanWind(obj)

  v = [obj.layer.windSpeed];
  out = sum( [obj.layer.fractionnalR0].*v.^(5./3) )^(3/5);

end % get.meanWind()

%% Get Greenwood frequency
function out = get.greenwoodFrequency(obj)

  v = [obj.layer.windSpeed];
  out = 0.4292*sum( [obj.layer.fractionnalR0].*v.^(5./3) )^(3/5)/obj.r0;

end % get.greenwoodFrequency()

function out = fourierPhaseScreen(atm,D,nPixel,nMap)
%% FOURIERPHASESCREEN Phase screen computation
%
% map = fourierPhaseScreen(atm,D,nPixel) Computes a square
% phase screen of D meter and sampled with nPixel using the
% Fourier method
%
% map = fourierPhaseScreen(atm) Computes a square phase screen
% of atm.layer.D meter and sampled with atm.layer.nPixel using
% the Fourier method; atm must contain only one turbulent
% layer!
%
% See also atmosphere

  if nargin<4
    nMap = 1;
  end
  if nargin<2
    D = atm.layer.D;
    nPixel = atm.layer.nPixel;
  end

  N = 4*nPixel;
  L = (N-1)*D/(nPixel-1);
  [fx,fy]  = freqspace(N,'meshgrid');
  [fo,fr]  = cart2pol(fx,fy);
  fr  = fftshift(fr.*(N-1)/L./2);
  clear fx fy fo
  psdRoot = sqrt(phaseStats.spectrum(fr,atm)); % Phase FT magnitude
  clear fr
  fourierSampling = 1./L;
  u = 1:nPixel;
  out = zeros(nPixel,nPixel,nMap);
  for kMap=1:nMap
    map = psdRoot.*fft2(randn(atm.rngStream,N))./N;
    % White noise filtering
    map = real(ifft2(map).*fourierSampling).*N.^2;
    out(:,:,kMap) = map(u,u);
  end

end % fourierPhaseScreen()

function map = choleskyPhaseScreen(atm,D,nPixel,nMap)
%% CHOLESKYPHASESCREEN Phase screen computation
%
% map = choleskyPhaseScreen(obj,D,nPixel) Computes a square
% phase screen of D meter and sampled with nPixel using the
% Cholesky factorisation of the atmospheric phase covariance
% matrix
%
% map = choleskyPhaseScreen(obj,D,nPixel,covarianceMatrix)
% Computes a square phase screen of D meter and sampled with
% nPixel using the Cholesky factorisation of the atmospheric
% phase covariance matrix covarianceMatrix
%
% See also chol and atmosphere

  if nargin<4
    nMap = 1;
  end
  if nargin==1
    for kLayer=1:atm.nLayer
      nPixel = atm.layer(kLayer).nPixel;
      if isempty(atm.layer(kLayer).choleskyFact)
        fprintf('Computing the Cholesky factor matrix!\n')
        D = atm.layer(kLayer).D;
        [x,y] = meshgrid(linspace(-1,1,nPixel)*D/2);
        atm.layer(kLayer).choleskyFact = chol( ...
        phaseStats.covarianceToeplitzMatrix(slab(atm,kLayer),complex(x,y)),...
                                            'lower');
      end
      map = atm.layer(kLayer).choleskyFact*randn(atm.rngStream,nPixel^2,nMap);
      map = reshape(map,nPixel,nPixel,nMap);
      atm.layer(kLayer).phase = map;
    end
  else
    fprintf('Computing the Cholesky factor matrix!\n')
    [x,y] = meshgrid(linspace(-1,1,nPixel)*D/2);
    choleskyFact = chol( ...
    phaseStats.covarianceToeplitzMatrix(atm,complex(x,y)) ,'lower');
    map = choleskyFact*randn(atm.rngStream,nPixel^2,nMap);
    map = reshape(map,nPixel,nPixel,nMap);
  end

end % choleskyPhaseScreen()


end % methods


methods (Static)


function out = r0VsWavelength(r0Wavelength,r0,wavelength)
%% R0VSWAVELENGTH r0 versus wavelength
%
% out = r0VsWavelength(r0Wavelength,r0,wavelength) computes the
% value of r0(r0Wavelength) at the new wavelength

  if ischar(r0Wavelength)
    r0Wavelength = eval(['photometry.',r0Wavelength]);
  end
  if ischar(wavelength)
    wavelength = eval(['photometry.',wavelength]);
  end
  out = r0.*(wavelength./r0Wavelength).^1.2;

end % r0VsWavelength()

end % methods (Static)


end % classdef