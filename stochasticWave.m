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
        nSample;
        samples;
        saveImage;
        % image 
        image
        % integratedImage
        cumImage=0;
        % phase variance buffer
        phaseVar;
        % logging flag
        logging = false;
    end
    
    properties (Dependent)
        % wave amplitude
        amplitude;
        % wave phase        
        phase;   
        % intensity
        intensity;
        % buffer sequence
        bufSeq
    end
    
    properties (Dependent,SetAccess=private)
        % the wave
        wave;
        % mean value removed phase map
        meanRmPhase; 
    end
    
    properties (Access=private)
        p_bufSeq;
        bufferLength = 10000;
        kBufSeq;
    end
    properties (Access=protected)
        p_amplitude = 1;
        p_phase = 0;
    end
    
    methods
        
        %% Constructor
        function obj = stochasticWave
        end
        
        %% Get/Set bufSeq
        function out = get.bufSeq(obj)
            out = obj.p_bufSeq;
        end
        function set.bufSeq(obj,val)
            obj.p_bufSeq = val;
            obj.kBufSeq = 1;
            obj.samples  = zeros(sum(obj.p_bufSeq),obj.bufferLength);
            obj.nSample  = 0;
        end
        
        %% Get/Set amplitude
        function out = get.amplitude(obj)
            out = obj.p_amplitude;
            out(~obj.mask) = 0;
        end
        function set.amplitude(obj,val)
            obj.p_amplitude = bsxfun( @times, obj.p_amplitude , val);
        end
    
        %% Get/Set phase
        function out = get.phase(obj)
            out = obj.p_phase;
            out(~obj.mask) = 0;
        end
        function set.phase(obj,val)
            obj.p_phase = bsxfun( @plus, obj.p_phase , val);
            if obj.logging
                buf = stats(obj,@var);
                if isempty(obj.phaseVar)
                    obj.phaseVar = zeros(obj.bufferLength,length(buf));
                    obj.kBufSeq = 0;
                end
                if length(buf)>size(obj.phaseVar,2)
                    obj.phaseVar = repmat( obj.phaseVar , 1 , length(buf));
                end
                obj.kBufSeq = obj.kBufSeq + 1;
                obj.phaseVar(obj.kBufSeq,:) = buf;
            end
%             if ~isempty(obj.p_bufSeq)
%                 if obj.p_bufSeq(obj.kBufSeq)
%                     obj.nSample = obj.nSample + 1;
%                     obj.samples(obj.nSample) = std(obj);
%                 end
%                 obj.kBufSeq = obj.kBufSeq + 1;
%                 if obj.kBufSeq>length(obj.p_bufSeq)
%                     obj.kBufSeq = 1;
%                 end
%             end
%             if obj.saveImage
%                 a = obj.wave;
%                 [n,m] = size(a);
%                 obj.image = abs( fft2( a , 2*n , 2*m ) ).^2;
%                 obj.cumImage = obj.cumImage + obj.image;
%             end
        end
        
        %% Get the wave
        function out = get.wave(obj)
            out = bsxfun( @times , obj.amplitude , exp(1i.*obj.phase) );
        end
        
        %% Get intensity
        function out = get.intensity(obj)
            out = obj.wave.*conj(obj.wave);
        end
        
        %% Get the meanRmPhase propertye
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
            %% RESET Reset wave properties
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
            %% CATAMPLITUDE Concatenate amplitudes
            %
            % out = catAmplitude(obj) concatenates the amplitudes of an
            % array of stochasticWave objects
            
            out = cat( ndims(obj) , obj.amplitude);
        end
        function out = catPhase(obj)
            %% CATPHASE Concatenate phases
            %
            % out = catPhase(obj) concatenates the phases of an array of
            % stochasticWave objects
            
            out = cat( ndims(obj) , obj.phase);
        end
        function out = phaseVector(obj)
            %% PHASEVECTOR Concatenate phases in 1 single vector
            %
            % out = phaseVector(obj) concatenates phase(s) values within
            % the mask in 1 single vector 
            
            nObj = length(obj);
            out = cell(nObj,1);
            for kSrc=1:nObj
%                 m_mask = repmat(obj(kSrc).mask,[1,1,size(obj(kSrc).phase,3)]);
                buf = utilities.toggleFrame(obj(kSrc).phase,2);
                out{kSrc} = buf(obj(kSrc).mask,:);
            end
            out = cell2mat(out);
        end
        function out = catMeanRmPhase(obj)
            %% CATPHASE Concatenate mean removed phases
            %
            % out = catMeanRmPhase(obj) concatenates the mean removed
            % phases of an array of stochasticWave objects
            
            out = cat( ndims(obj) , obj.meanRmPhase);
        end
        function out = catWave(obj)
            %% CATWAVE Concatenate waves
            %
            % out = catWave(obj) concatenates the waves of an array of
            % stochasticWave objects
            out = cat( ndims(obj) , obj.wave);
        end
        
        function out = mean(obj)
            %% MEAN Average or mean value
            %
            % out = mean(obj) computes the average or mean value of the
            % amplitude and/or the phase within the mask. If either
            % amplitude or phase is a stochastic property out is a single
            % value; if both amplitude and phase are stochastic out is a 2
            % element vector with amplitude and phase results in this order

            out = stats(obj,@mean);
        end
        
        function out = var(obj)
            %% VAR Variance
            %
            % out = var(obj) computes the variance of the amplitude and/or
            % the phase within the mask. If either amplitude or phase is a
            % stochastic property out is a single value; if both amplitude
            % and phase are stochastic out is a 2 element vector with
            % amplitude and phase results in this order

            out = stats(obj,@var);
       end
        
        function out = std(obj)
            %% STD Standard deviation
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
%                 k = 1;
%                 if obj.stochasticAmplitude
%                     out(k) = fun(obj.amplitude(obj.mask));
%                     k = k+1;
%                 end
                if obj.stochasticPhase
                    buf = utilities.toggleFrame(obj.phase,2);
                    out = fun(buf(obj.mask,:));
                end
            end
        end
        
    end
    
end