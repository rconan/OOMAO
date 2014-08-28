classdef abcd < handle
    % Create a abcd object
    
    properties
        % ABCD matrix
        matrix;    
        % input offset and angle ray vector
        offsetAngle;
        % element thickness
        thickness = 0;
        % element offset
        offset;
        % element stop width
        stopWidth;
        % element stop offset
        stopOffset;
        % z-axis propagation direction (1: forward, -1:backward)
        zPropDir = 1;
    end
    
    properties (Abstract)
        % tag
        tag;
    end        

    methods
        
        %% Constructor
        function obj = abcd(varargin)
            
            p = inputParser;
            p.addParameter('offset',0,@isnumeric);
            p.addParameter('stopWidth',Inf,@isnumeric);
            p.addParameter('stopOffset',Inf,@isnumeric);
            p.addParameter('zPropDir',1,@isnumeric);
            p.parse(varargin{:});

            obj.offset = p.Results.offset;
            obj.stopWidth = p.Results.stopWidth;
            obj.stopOffset = p.Results.stopOffset;
            obj.zPropDir = p.Results.zPropDir;
            
        end
        
        function relay(obj,src)
            %% RELAY source propagation through element
            
            % input offset angle vector
            src.offsetAngle(1,:) = src.offsetAngle(1,:)  - obj.offset;
            obj.offsetAngle = src.offsetAngle;
            % update abcd matrix
            src.offsetAngle = obj.matrix*src.offsetAngle;
        end
    end
    
end