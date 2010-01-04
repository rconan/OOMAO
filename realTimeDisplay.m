%% realTimeDisplay Class Definition
% The realTimeDisplay class refreshes the graphic display of other objects
% at a user defined rate. A realTimeDisplay display is created with:
% obj = realTimeDisplay(callback) where callback is the handle of the
% graphic display function. 
% For example, to display a CCD frame at the default rate of 1Hz, ypu can
% do the following:
% rtd = realTimeDisplay( @(x) imagesc(ccd) );
% start(rtd.paceMaker);
classdef realTimeDisplay < handle
%%

    %% Properties
    properties
        displayRate; % Hz (default: 1Hz)
        paceMaker;
        callBack;
    end
    
    %% Methods
    methods
        
        %%
        % * Constructor
        function obj = realTimeDisplay(callback)
            error(nargchk(1, 1, nargin))
            obj.callBack = callback;
            obj.paceMaker = timer;
            obj.paceMaker.name = 'Real Time Display';
            obj.paceMaker.TimerFcn = @timerCallBack;
            obj.paceMaker.ExecutionMode = 'FixedDelay';
            function timerCallBack( timerObj, event)
                obj.callBack();
            end
            obj.displayRate = 1;
        end
        
        %%
        % * Destructor
        function delete(obj)
            if isvalid(obj.paceMaker)
                if strcmp(obj.paceMaker.Running,'on')
                    stop(obj.paceMaker)
                end
                delete(obj.paceMaker)
            end
%             logBook.add(obj,'Terminated!')
        end

        %%
        % * Interactive display rate setting
        function obj = set.displayRate(obj, val)
            obj.displayRate = val;
            if strcmp(obj.paceMaker.Running,'on')
                stop(obj.paceMaker)
                obj.paceMaker.Period = 1/val;
                start(obj.paceMaker)
            else
                obj.paceMaker.Period = 1/val;
            end
        end


    end

end