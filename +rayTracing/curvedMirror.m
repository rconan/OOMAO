classdef curvedMirror < rayTracing.abcd
    % Create a curvedMirror object
    %
    % tl = curvedMirror(focalLenght) creates a curvedMirror object from the
    % radius of curvature
    %
    % tl = curvedMirror(...,'offset',offset)
    
    properties
        % thin lens parameters vector
        params;
        % tag
        tag = 'freeSpace';
    end
    
    methods
        
        %% Constructor
        function obj = curvedMirror(params,varargin)
            
            obj = obj@rayTracing.abcd(varargin{:});
            
            obj.params = params;
            obj.matrix = [1 0 ; -2/params 1];
            obj.zPropDir = -1;            
            
        end
        
        function relay(obj,src)
            relay@rayTracing.abcd(obj,src);
            obj.thickness = sqrt(obj.params.^2 - obj.offsetAngle(1,:).^2) - obj.params;
        end
        
    end
    
end