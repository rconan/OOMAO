classdef thinLens < rayTracing.abcd
    % Create a thinLens object
    %
    % tl = thinLens(focalLenght) creates a thinLens object from the lens
    % focal length
    %
    % tl = thinLens(...,'offset',offset)
    
    properties
        % thin lens parameters vector
        params;
        % tag
        tag = 'freeSpace';
    end
    
    methods
        
        %% Constructor
        function obj = thinLens(params,varargin)
            
            obj = obj@rayTracing.abcd(varargin{:});
            
            obj.params = params;
            obj.matrix = [1 0 ; -1/params 1];
            
        end
        
    end
    
end