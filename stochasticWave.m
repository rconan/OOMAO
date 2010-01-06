classdef stochasticWave < handle
    % STOCHASTICWAVE Create a stochasticWave object

    properties
        % mask to retain amplitude and phase values used for statistics
        % moments evaluation
        mask;
        % true if the phase is a stochastic property
        stochasticPhase = true;
        % true if the amplitude is a stochastic property
        stochasticAmplitude = false;
    end
    
    properties (Dependent)
        % wave amplitude
        amplitude;
        % wave phase        
        phase;        
    end
    
    properties (Dependent,SetAccess=private)
        % the wave
        wave;
        % mean value removed phase map
        meanRmPhase; 
    end
    
    properties (Access=private)
        p_amplitude = 1;
        p_phase = 0;
    end
    
    methods
        
        % Constructor
        function obj = stochasticWave
        end
        
        % Get/Set amplitude
        function out = get.amplitude(obj)
            out = obj.p_amplitude;
            out(~obj.mask) = 0;
        end
        function set.amplitude(obj,val)
            obj.p_amplitude = bsxfun( @times, obj.p_amplitude , val);
        end
        % Get/Set phase
        function out = get.phase(obj)
            out = obj.p_phase;
            out(~obj.mask) = 0;
        end
        function set.phase(obj,val)
            obj.p_phase = bsxfun( @plus, obj.p_phase , val);
        end
        
        % Get the wave
        function out = get.wave(obj)
            out = bsxfun( @times , obj.amplitude , exp(1i.*obj.phase) );
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
        
        function varargout = reset(obj)
            % RESET Reset wave properties
            %
            % reset(obj) resets the mask to [], the amplitude to 1 and the
            % phase to 0
            %
            % obj = reset(obj) resets and returns the resetted object
            
            for kObj = 1:numel(obj);
                obj(kObj).mask        = [];
                obj(kObj).p_amplitude = 1;
                obj(kObj).p_phase     = 0;
            end
            if nargout>0
                varargout{1} = obj;
            end
        end
        
        function out = catAmplitude(obj)
            % CATAMPLITUDE Concatenate amplitudes
            %
            % out = catAmplitude(obj) concatenates the amplitudes of an
            % array of stochasticWave objects
            
            out = cat( ndims(obj) , obj.amplitude);
        end
        function out = catPhase(obj)
            % CATPHASE Concatenate phases
            %
            % out = catPhase(obj) concatenates the phases of an array of
            % stochasticWave objects
            
            out = cat( ndims(obj) , obj.phase);
        end
        function out = catWave(obj)
            % CATWAVE Concatenate waves
            %
            % out = catWave(obj) concatenates the waves of an array of
            % stochasticWave objects
            out = cat( ndims(obj) , obj.wave);
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