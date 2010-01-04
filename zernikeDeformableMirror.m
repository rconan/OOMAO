classdef zernikeDeformableMirror < telescope & zernike
    
    properties
        height;
        conjD;
    end
    
    methods
        
        % Constructor
        function obj = zernikeDeformableMirror(zernParam,telParam,height)
            obj = obj@telescope(telParam{:});
            obj = obj@zernike(zernParam{:});
            obj.height = 0;
            if nargin>2
                obj.height = height;
            end
            obj.conjD = obj.diameterAt(obj.height);
        end
        
        function varargout = footprintProjection(obj,src)
            delta = obj.height.*tan(src.zenith).*...
                [cos(src.azimuth),sin(src.azimuth)];
            delta = delta*2/obj.conjD;
            alpha = obj.conjD./obj.D;
            P = obj.smallFootprintExpansion(delta,alpha);
            varargout{1} = P;
            if nargout>1
                o = linspace(0,2*pi,101);
                varargout{2} = cos(o)./alpha + delta(1);
                varargout{3} = sin(o)./alpha + delta(2);
            end
        end
        
    end
    
end