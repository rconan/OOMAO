classdef constants
   
    properties (Constant=true)
        radian2arcsec = 180*3600/pi;
        radian2mas    = 1e3*180*3600/pi;
        radian2arcmin = 180*60/pi;
        arcsec2radian = pi/180/3600;
        arcmin2radian = pi/180/60;
        plank         = 6.62606896e-34; % J.s=kg.m^2/s
        c             = 299792458; % light of speed m/s
        Me            = 5.9721986e24; % earth mass [kg]
        Re            = 6378.14e3; % earth equatorial radius [m]
        G             = 6.67e-11; % Newtonian gravitational constant [Nm^2/kg^2] (1N=1kg.m.s^-2)
    end
    
    methods (Static)
        
        function radian = arcmin(val)
            radian = val*constants.arcmin2radian;
        end
        
    end
    
end