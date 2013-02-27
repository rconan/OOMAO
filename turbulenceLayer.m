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
        sampling;
        zern;
        choleskyFact;
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
        
        function host2gpu(obj)
            nLayer = length(obj);
            for kLayer=1:nLayer
                obj(kLayer).altitude       = gsingle( obj(kLayer).altitude );
                obj(kLayer).fractionnalR0  = gsingle( obj(kLayer).fractionnalR0 );
                obj(kLayer).windSpeed      = gsingle( obj(kLayer).windSpeed );
                obj(kLayer).windDirection  = gsingle( obj(kLayer).windDirection );
                obj(kLayer).amplitude      = gsingle( obj(kLayer).amplitude );
                obj(kLayer).D              = gsingle( obj(kLayer).D );
                obj(kLayer).nPixel         = gsingle( obj(kLayer).nPixel );
                obj(kLayer).sampling       = gsingle( obj(kLayer).sampling );
                obj(kLayer).phase          = gsingle( obj(kLayer).phase );
            end
        end
        
        function gpu2host(obj)
            nLayer = length(obj);
            for kLayer=1:nLayer
                obj(kLayer).altitude       = double( obj(kLayer).altitude );
                obj(kLayer).fractionnalR0  = double( obj(kLayer).fractionnalR0 );
                obj(kLayer).windSpeed      = double( obj(kLayer).windSpeed );
                obj(kLayer).windDirection  = double( obj(kLayer).windDirection );
                obj(kLayer).amplitude      = double( obj(kLayer).amplitude );
                obj(kLayer).D              = double( obj(kLayer).D );
                obj(kLayer).nPixel         = double( obj(kLayer).nPixel );
                obj(kLayer).sampling       = double( obj(kLayer).sampling );
                obj(kLayer).phase          = double( obj(kLayer).phase );
            end
        end
        
    end
    
end