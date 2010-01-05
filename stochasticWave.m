classdef stochasticWave < handle
    % STOCHASTIC WAVE Create a stochasticWave object

    properties
        % mask to retain amplitude and phase values used for statistics
        % moments evaluation
        mask;
        % wave amplitude
        amplitude;
        % wave phase        
        phase;
        % true if the phase is a stochastic property
        stochasticPhase = true;
        % true if the amplitude is a stochastic property
        stochasticAmplitude = false;
    end
    
    properties (Dependent,SetAccess=private)
        % the wave
        wave;
        % mean value removed phase map
        meanRmPhase; 
    end

    
    methods
        
        % Constructor
        function obj = stochasticWave
        end
        
        % Get the wave
        function out = get.wave(obj)
            out = obj.amplitude.*exp(1i.*obj.phase);
        end
        
        % Get the meanRmPhase propertye
        function meanRmPhase = get.meanRmPhase(obj)
            buffer = obj.stochasticAmplitude;
            obj.stochasticAmplitude = false; % Set to false in order for the mean method to return only phase value
            meanPhase   = mean(obj); % Set the mask if empty
            meanRmPhase = obj.phase;
            meanRmPhase(obj.mask) ...
                        = obj.phase(obj.mask) - meanPhase;
            obj.stochasticAmplitude = buffer;
        end
        
        function out = mean(obj)
            % MEAN Average or mean value
            %
            % out = mean(obj) computes the average or mean value of the
            % amplitude and/or the phase within the mask. If either
            % amplitude or phase is a stochastic property out is a single
            % value; if both amplitude and phase are stochastic out is a 2
            % element vector with amplitude and phase results in this order

            out = stats(obj,@mean);
        end
        
        function out = var(obj)
            % VAR Variance
            %
            % out = var(obj) computes the variance of the amplitude and/or
            % the phase within the mask. If either amplitude or phase is a
            % stochastic property out is a single value; if both amplitude
            % and phase are stochastic out is a 2 element vector with
            % amplitude and phase results in this order

            out = stats(obj,@var);
       end
        
        function out = std(obj)
            % STD Standard deviation
            %
            % out = std(obj) computes the standard deviation of the
            % amplitude and/or the phase within the mask. If either
            % amplitude or phase is a stochastic property out is a single
            % value; if both amplitude and phase are stochastic out is a 2
            % element vector with amplitude and phase results in this order
            
            out = stats(obj,@std);
        end
        
    end
        
    methods (Access=private)
        
        function out = stats(obj,fun)
            if numel(obj)>1
                out = cellfun( @(x) stats(x,fun), num2cell(obj) );
            else
                if isempty(obj.mask)
                    disp(' @(stochasticWave)> Setting the mask!')
                    obj.mask = ones(size(obj.amplitude));
                end
                k = 1;
                if obj.stochasticAmplitude
                    out(k) = fun(obj.amplitude(obj.mask));
                    k = k+1;
                end
                if obj.stochasticPhase
                    out(k) = fun(obj.phase(obj.mask));
                end
            end
        end
        
    end
    
end