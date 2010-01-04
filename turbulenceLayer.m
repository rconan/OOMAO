% TURBULENCELAYER Create a turbulence layer object array
%
% turb = turbulenceLayer(altitude,fractionnalR0) creates a turbulence layer
% object array from the atitudes and the fractionnalR0 of the layers
%
% turb = turbulenceLayer(altitude,fractionnalR0,windSpeed,windDirection)
% creates a turbulence layer object array from the atitudes, the
% fractionnalR0, the wind speeds and the wind directions of the layers

classdef turbulenceLayer < handle
    
    properties
        altitude;
        fractionnalR0;
        windSpeed;
        windDirection;
        amplitude;
        D;
        nPixel;
        zern;
    end
    
    properties (Dependent=true)
        wave;
    end
        
    properties (SetObservable=true)
        % the image at each lenslet focus
        phase;
    end
   
    methods
        
        % Constructor
        function obj = turbulenceLayer(...
                altitude,...
                fractionnalR0,...
                windSpeed,...
                windDirection)
            if nargin~=0
                error(nargchk(1, 4, nargin))
                nLayer = length(altitude);
                obj(nLayer) = turbulenceLayer;
                for kLayer=1:nLayer
                    obj(kLayer).altitude      = altitude(kLayer);
                    obj(kLayer).fractionnalR0 = fractionnalR0(kLayer);
                    if nargin>2
                        obj(kLayer).windSpeed     = windSpeed(kLayer);
                        obj(kLayer).windDirection = windDirection(kLayer);
                    end
                end
            end
        end
        
    end
    
end