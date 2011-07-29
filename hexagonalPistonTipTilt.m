classdef hexagonalPistonTipTilt < handle
    % INFLUENCEFUNCTION Create an influence function object
    %
    % obj = influenceFunction('monotonic',mech) creates a cubic Bezier
    % influence function monotically decreasing from 1 to 0 with the
    % mechanical coupling mech
    %
    % obj = influenceFunction('overshoot',mech) creates a cubic Bezier
    % influence function with a negative overshoot and the mechanical
    % coupling mech
    %
    % Try show(obj) to view a cut of the influence function
    
    properties
        % modes
        modes
        % influence function tag
        tag = 'HEXAGONAL PISTON TIP TILT MODES';
    end
     
    properties (Access=private)
        log
    end
    
    methods
        
        %% Constructor
        function obj = hexagonalPistonTipTilt()
            obj.log = logBook.checkIn(obj);
        end
                
        %% Destructor
        function delete(obj)
            checkOut(obj.log,obj)
        end
        
        function out = mtimes(obj,c)
            %% MTIMES Influence function multiplication
            %
            % v = obj*c multiply the influence functions in obj by the
            % vector c
            
            out = obj.modes*c;
        end

        function setInfluenceFunction(obj,nIF,resolution,validActuator,~,~)
            nIF = nIF/3;
            nCycle = roots([3,3,1-nIF]);
            nCycle(nCycle<0) = [];
            [~,center] = utilities.hexagonalArray(nCycle);
            u0 = linspace(-1,1,resolution)*(nCycle-1);
            pitch = resolution/(2*(nCycle-1));
            [tip,tilt] = meshgrid(u0*pitch);
            tip  = tip(:);
            tilt = tilt(:);
            center = center*pitch;
            nValid = sum(validActuator(:));
            obj.modes = zeros(resolution^2,nValid*3);
            count = 0;
            for kValid=1:nIF
                if validActuator(kValid)
                    xc = real(center(kValid));
                    yc = imag(center(kValid));
                    buf = utilities.piston(pitch*2/sqrt(3),resolution,...
                        xc,yc,'shape','hex');
                    buf = buf(:);
                    count = count + 1;
                    obj.modes(:,count) = buf;
                    count = count + 1;
                    obj.modes(:,count) = 2*buf.*(tip-xc)/pitch;
                    count = count + 1;
                    obj.modes(:,count) = 2*buf.*(tilt-yc)/pitch;
                end
            end
        end
        
    end
    
end