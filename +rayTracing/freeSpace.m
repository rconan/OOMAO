classdef freeSpace < rayTracing.abcd
    % Create a freeSpace object
    %
    % tl = freeSpace(focalLenght) creates a freeSpace object from the lens
    % focal length
    %
    % tl = freeSpace(...,'offset',offset)
    
    properties
        % thin lens parameters vector
        params;
        % tag
        tag = 'freeSpace';
    end
    
    methods
        
        %% Constructor
        function obj = freeSpace(params,varargin)
            
            obj = obj@rayTracing.abcd(varargin{:});
            
            obj.params = params;
            obj.matrix = [1 params ; 0 1];
            obj.thickness = params;
            
        end
        
    end
    
end