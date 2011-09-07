classdef skyAngle < double
    
    properties 
       tag = 'skyAngle'
       % angle
       angle
       % unit
       unit
    end
    
    methods
        
        %% Constructor
        function obj = skyAngle(angle_,unit_)
            
            if nargin<2
                unit_ = 'radian';
            end
            
            switch unit_
                case 'arcmin'
                    angle_ = angle_/constants.radian2arcmin;
                case 'arcsec'
                    angle_ = angle_/constants.radian2arcsec;
                case 'mas'
                    angle_ = angle_/constants.radian2mas;
                case 'degree'
                    angle_ = angle_*pi/180;
            end
            
            obj = obj@double(angle_);
            obj.angle = angle_;
            obj.unit  = unit_;
            
        end
        
        function display(obj)
            fprintf('  sky angle: %g %s\n',convert(obj,obj.unit),obj.unit)
        end
        
        function c = plus(a,b)
            c = skyAngle(a.angle + b.angle());
        end
        
        function out = RADIAN(obj)
            out = convert(obj,'radian');
        end
        
        function out = ARCMIN(obj)
            out = convert(obj,'arcmin');
        end
        
        function out = ARCSEC(obj)
            out = convert(obj,'arcsec');
        end
        
        function out = MAS(obj)
            out = convert(obj,'mas');
        end
        
        function out = DEGREE(obj)
            out = convert(obj,'degree');
        end
        
        function out = convert(obj,newUnit)
            switch newUnit
                case 'arcmin'
                    out = obj.angle*constants.radian2arcmin;
                case 'arcsec'
                    out = obj.angle*constants.radian2arcsec;
                case 'mas'
                    out = obj.angle*constants.radian2mas;
                case 'degree'
                    out = obj.angle*180/pi;
                otherwise
                    out = obj.angle;
            end
            obj.unit = newUnit;
        end
        
    end
    
end