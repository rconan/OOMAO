classdef controller < handle
% CONTROLLER Create a temporal controller object
%
% ctrlr = controller(wfs,dm,commandMatrix,nIteration) creates a
% temporal controller object from a wavefront sensor object, a
% deformable mirror object, the commandMatrix and the number of
% iterations. The controller is a closed loop integrator with no delay.
% The integrator gain is 0.5.
%
% ctrlr = controller(...,'tipTiltCompensator',ttMirror,'tipTiltCommandMatrix',ttM)
% Same as above but adding a separated tip-tilt mirror. The controller
% is a closed loop integrator with no delay. The integrator gain is
% 0.5.

% ctrlr = controller(...,'loopGain',gain) set the closed-loop
% integrator gain

properties
  % sensor object
  sensor;
  sensorData;
  sensorListener;
  nSensor;
  % tip-tilt sensor object
  tipTiltSensor;
  tipTiltSensorData;
  tipTiltSensorListener;
  nTipTiltSensor;
  % compensator object
  compensator;
  compensatorData;
  nCompensator;
  % tip=tilt compensator object
  tipTiltCompensator;
  tipTiltCompensatorData;
  nTipTiltCompensator;
  tipTiltDelay;
  % input phase
  input;
  % input variance
  inputVar;
  % output phase
  output;
  % output variance
  outputVar;
  % closed-loop or open-loop
  type;
  % delay between sensing and compensating
  delay
  % command matrix
  M
  % tip-tilt command matrix
  tipTiltM
  % integrator gain
  gain;
  % controller temporal steps
  nIteration;
  kIteration=0;
  % calibration matrix D and SVD
  D;
  U;
  S;
  V;
  % tag
  tag = 'CONTROLLER';
end % properties


properties (Access=private)
  log;
end % properties (Access=private)


methods


%% Constructor
function obj = controller(sensor,compensator,commandMatrix,nIteration,varargin)

  % Parsing input parameters
  p = inputParser;
  p.addRequired('sensor', @(x) isa(x,'shackHartmann') );
  p.addRequired('compensator', @(x) isa(x,'deformableMirror') );
  p.addRequired('nIteration', @isnumeric );
  p.addOptional('commandMatrix', [] ,@isnumeric );
  p.addParamValue('nSensor', 1, @isnumeric );
  p.addParamValue('nTipTiltSensor', 1, @isnumeric );
  p.addParamValue('nCompensator', 1, @isnumeric );
  p.addParamValue('type', 'closedLoop', @ischar );
  p.addParamValue('delay', 0, @isnumeric );
  p.addParamValue('gain', 0.5, @isnumeric );
  p.addParamValue('tipTiltSensor', [] , @(x) isa(x,'shackHartmann') );
  p.addParamValue('tipTiltCompensator', [] , @(x) isa(x,'deformableMirror') );
  p.addParamValue('tipTiltCommandMatrix', [] , @isnumeric );
  p.addParamValue('tipTiltDelay', 0 , @isnumeric );
  p.parse(sensor,compensator,commandMatrix,nIteration,varargin{:});

  % Allocating properties
  obj.sensor        = p.Results.sensor;
  obj.tipTiltSensor = p.Results.tipTiltSensor;
  obj.nSensor       = p.Results.nSensor;
  obj.nTipTiltSensor ...
                    = p.Results.nTipTiltSensor;
  obj.compensator   = p.Results.compensator;
  obj.tipTiltCompensator ...
                    = p.Results.tipTiltCompensator;
  obj.nCompensator  = p.Results.nCompensator;
  obj.M             = p.Results.commandMatrix;
  obj.tipTiltM      = p.Results.tipTiltCommandMatrix;
  obj.nIteration    = p.Results.nIteration;
  obj.type          = p.Results.type;
  obj.delay         = p.Results.delay;
  obj.gain          = p.Results.gain;
  obj.tipTiltDelay  = p.Results.tipTiltDelay;

  nSrc = 1;
  obj.inputVar      = zeros(obj.nIteration,nSrc);
  obj.outputVar     = zeros(obj.nIteration,nSrc);
  obj.sensorData    = ...
  zeros(obj.sensor.nSlope,obj.nSensor,obj.nIteration+obj.delay+1);
  obj.compensatorData = ...
  zeros(obj.compensator.nValidActuator,obj.nIteration+obj.delay+1);
  obj.tipTiltSensorData = ...
  zeros(2,obj.nTipTiltSensor,obj.nIteration+obj.delay+1);
  obj.tipTiltCompensatorData = ...
  zeros(2,obj.nTipTiltSensor,obj.nIteration+obj.delay+1);

  obj.log = logBook.checkIn(obj);

  % Listener setting
  if isempty(obj.tipTiltCompensator)
    obj.sensorListener = ...
    addlistener(obj.sensor,'slopes','PostSet',@obj.closedLoop);
    add(obj.log,obj,'Closed-loop integrator with:\n . 1 DM,\n . 1 WFS!')
  elseif isempty(obj.tipTiltSensor)
    obj.sensorListener = ...
    addlistener(obj.sensor,'slopes','PostSet',@obj.closedLoopWithtTipTiltComp);
    add(obj.log,obj, ...
    'Closed-loop integrator with:\n . 1 DM,\n . 1 TT mirror\n . 1 WFS!')
  else
    obj.sensorListener = addlistener(obj.sensor,'slopes','PostSet',...
                                     @obj.closedLoopWithtTipTiltSensComp);
    add(obj.log,obj,...
    'Closed-loop integrator with:\n . 1 DM,\n . 1 TT mirror,\n . 1 WFS,\n . 1 TT WFS!')
  end
  obj.sensorListener.Enabled = false;

end % controller()


% Destructor
function delete(obj)

  if ~isempty(obj.log)
    checkOut(obj.log,obj)
  end

end % delete()


function calibration(obj,gs,nTrunc)
%% CALIBRATION DM to WFS calibration
%
% calibration(ctrlr,gs.*tel)
% calibration(ctrlr,gs.*tel,nTrunc)

  if isempty(obj.D)
    add(obj.log,obj,'Computing the poke matrix')

    dm = obj.compensator;
    wfs = obj.sensor;

    fprintf(' ___ CALIBRATION ___\n')
    calibDmCommands = speye(dm.nValidActuator)*gs.wavelength/4;
    if dm.nValidActuator>1000
      steps           = 40;
    else
      steps = 1;
    end

    nC              = floor(dm.nValidActuator/steps);
    u               = 0;
    obj.D           = zeros(wfs.nSlope,dm.nValidActuator);
    gs  = gs*dm*wfs;
    buf = dm.coefs;

    fprintf(' . actuators range:          ')
    while u(end)<dm.nValidActuator
      u = u(end)+1:min(u(end)+nC,dm.nValidActuator);
      fprintf('\b\b\b\b\b\b\b\b\b%4d:%4d',u(1),u(end))
      dm.coefs = calibDmCommands(:,u);
      +gs;
      obj.D(:,u) = wfs.slopes/(gs.wavelength/4);
    end
    fprintf('\n--------------------\n')


%                 dm.coefs = eye(dm.nValidActuator)*gs.wavelength/4;
%                 gs = gs*dm*wfs;
%                 obj.D = wfs.slopes/(gs.wavelength/4);

    dm.coefs = buf;
    gs = gs.*gs.opticalPath{1};

    figure
    subplot(1,2,1)
    imagesc(obj.D)
    axis equal tight
    colorbar

    [obj.U,obj.S,obj.V] = svd(obj.D);
    s = diag(obj.S);

    subplot(1,2,2)
    semilogy(s,'.')
    drawnow

  end

  if isempty(obj.tipTiltM)

    if ~isempty(obj.tipTiltCompensator)

      tt = obj.tipTiltCompensator;
      buf = tt.coefs;
      tt.coefs = eye(2)*gs.wavelength/4;

      if isempty(obj.tipTiltSensor)

        wfs = obj.sensor;
        buf1 = wfs.rmMeanSlopes;
        wfs.rmMeanSlopes = true;
        gs = gs*tt*wfs;
        Dtt = wfs.meanSlopes/(gs.wavelength/4);
        obj.tipTiltM = pinv(Dtt);
        wfs.rmMeanSlopes = buf1;

      else

        ttWfs = obj.tipTiltSensor;
        gs = gs*tt*ttWfs;
        DttWfs = ttWfs.slopes/(gs.wavelength/4);
        obj.tipTiltM = pinv(DttWfs);

      end
      tt.coefs = buf;
      gs = gs.*gs.opticalPath{1};

    end

  end

  %             nTrunc = 4;
  if nargin>2
    s = diag(obj.S);
    %index = s/s(1)>=nTrunc;%length(s)-nTrunc;
    %obj.M = obj.V(:,1:index)*diag(1./s(1:index))*obj.U(:,1:index)';
    index = s/s(1)>nTrunc;%length(s)-nTrunc;
    fprintf(' Number of thresholded values: %d out of %d\n',...
            sum(~index),length(index))
    obj.M = obj.V(:,index)*diag(1./s(index))*obj.U(:,index)';
  end

end % calibration()


%function relay(obj,src)

%%             nSrc = length(src);
     %obj.kIteration = obj.kIteration + 1;

     %if strcmpi(obj.type,'closedLoop') % CLOSED-LOOP

         %obj.inputVar(obj.kIteration) = var(src);
         %src = src*obj.compensator*obj.sensor;
         %obj.outputVar(obj.kIteration) = var(src);

         %obj.sensorData(:,:,obj.kIteration+obj.delay+1)      = obj.sensor.slopes;
         %obj.compensatorData(:,obj.kIteration+obj.delay+1) = obj.compensator.coefs;

         %obj.compensator.coefs = obj.compensator.coefs - obj.gain*obj.M*obj.sensorData(:,:,obj.kIteration);

     %elseif strcmpi(obj.type,'openLoop') % OPEN-LOOP

         %obj.inputVar(obj.kIteration,:) = var(src);
         %src = src*obj.sensor*obj.compensator;
         %obj.outputVar(obj.kIteration,:) = var(src);

         %obj.sensorData(:,:,obj.kIteration+obj.delay+1)      = obj.sensor.slopes;
         %obj.compensatorData(:,:,obj.kIteration+obj.delay+1) = obj.compensator.coefs;

         %obj.compensator.coefs = obj.M*obj.compensatorData(:,:,obj.kIteration);

     %else % OUPS!

         %warning('oomao:controller:relay','The controller type must be either closedLoop or openLoop!')

     %end

%end


function run(obj,runFun)

  switch class(runFun)
    case 'cell'
      for m_kIteration=1:obj.nIteration
        cellfun(@(f) f(), runFun)
      end
    case 'functionHandle'
      for m_kIteration=1:obj.nIteration
        runFun()
      end
    otherwise
      error('oomao:controller:run','run argument is either one function handle of a cell array of function handle!')
  end

end % run()


function closedLoop(obj,~,~) % obj, src, event
%% CLOSED LOOP closed loop integrator

  obj.kIteration = obj.kIteration + 1;
  obj.sensorData(:,:,obj.kIteration+obj.delay+1)    = obj.sensor.slopes;
  obj.compensatorData(:,obj.kIteration+obj.delay+1) = obj.compensator.coefs;

  obj.compensator.coefs = obj.compensator.coefs - ...
                          obj.gain*obj.M*obj.sensorData(:,:,obj.kIteration);

end % closedLoop()


function closedLoopWithtTipTiltComp(obj,~,~) % obj, src, event
%% CLOSEDLOOPWITHTTIPTILTCOMP closed loop integrator

  obj.kIteration = obj.kIteration + 1;

  % WFS
  obj.sensorData(:,:,obj.kIteration+obj.delay+1)           = ...
  obj.sensor.slopes;
  % TT WFS
  obj.tipTiltSensorData(:,:,obj.kIteration+obj.delay+1)    = ...
  obj.sensor.meanSlopes;
  % DM
  obj.compensatorData(:,obj.kIteration+obj.delay+1)        = ...
  obj.compensator.coefs;
  % TT Mirror
  obj.tipTiltCompensatorData(:,obj.kIteration+obj.delay+1) = ...
  obj.tipTiltCompensator.coefs;

  obj.compensator.coefs                                    = ...
  obj.compensator.coefs - obj.gain*obj.M*obj.sensorData(:,:,obj.kIteration);

  obj.tipTiltCompensator.coefs                             = ...
  obj.tipTiltCompensator.coefs - ...
  obj.gain*obj.tipTiltM*obj.tipTiltSensorData(:,:,obj.kIteration);

end % closedLoopWithtTipTiltComp()


function closedLoopWithtTipTiltSensComp(obj,~,~) % obj, src, event
%% CLOSEDLOOPWITHTTIPTILTSENSCOMP closed loop integrator

  obj.kIteration = obj.kIteration + 1;

  % WFS
  obj.sensorData(:,:,obj.kIteration+obj.delay+1) = obj.sensor.slopes;
  % TT WFS
  obj.tipTiltSensorData(:,:,obj.kIteration+obj.delay+1) = ...
  obj.tipTiltSensor.slopes;

  obj.compensator.coefs = ...
  obj.compensator.coefs - obj.gain*obj.M*obj.sensorData(:,:,obj.kIteration);

  if obj.kIteration>=obj.tipTiltDelay
    obj.tipTiltCompensator.coefs = ...
    bsxfun( @minus , obj.tipTiltCompensator.coefs , ...
    obj.gain*obj.tipTiltM*obj.tipTiltSensorData(:,:,obj.kIteration) );
  end

%% TODO: Fix below
%% DM
% obj.compensatorData(:,obj.kIteration+obj.delay+1) = obj.compensator.coefs;
%% TT Mirror
% obj.tipTiltCompensatorData(:,:,obj.kIteration+obj.delay+1) =...
%%obj.tipTiltCompensator.coefs;

end % closedLoopWithtTipTiltSensComp()


function openLoop(obj,~,~) % obj, src, event

  obj.kIteration = obj.kIteration + 1;
  obj.sensorData(:,:,obj.kIteration+obj.delay+1)    = obj.sensor.slopes;
  obj.compensatorData(:,obj.kIteration+obj.delay+1) = obj.compensator.coefs;
  obj.compensator.coefs = obj.M*obj.compensatorData(:,:,obj.kIteration);

end % openLoop()


end % methods


end % classdef