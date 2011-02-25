classdef photometry < handle
    
    properties
        wavelength % [micron]
        bandwidth  % [micron]
        zeroPoint  % [ph/s]
    end
    
    properties (Dependent)
        nPhoton   % # of photons
        magnitude % magnitude
    end
    
    properties (Access=private)
        p_nPhoton
        p_magnitude
    end
    
    methods
        
        %% Constructor
        function obj = photometry(w,bw,zp)
            obj.wavelength = w;
            obj.bandwidth  = bw;
            obj.zeroPoint  = zp/368;
        end
        
        %% Get and Set magnitude
        function out = get.magnitude(obj)
            out = obj.p_magnitude;
        end
        function set.magnitude(obj,val)
            obj.p_magnitude = val;
            obj.p_nPhoton = obj.zeroPoint*10^(-0.4*val);
%             fprintf(' @(photometry)> # of photon s^{-1}: %g\n',obj.p_nPhoton)
        end
        
        %% Get and Set nPhoton
        function out = get.nPhoton(obj)
            out = obj.p_nPhoton;
        end
        function set.nPhoton(obj,val)
            obj.p_nPhoton = val;
            obj.p_magnitude = -2.5*log10(val/obj.zeroPoint);
            fprintf(' @(photometry)> magnitude %4.2f\n',obj.p_magnitude)
        end
        
        function obj = plus(obj1,obj2)
            w   = obj1.zeroPoint*obj1.wavelength + obj2.zeroPoint*obj2.wavelength;
            bw  = obj1.bandwidth + obj2.bandwidth;
            zp  = obj1.zeroPoint + obj2.zeroPoint;
            obj = gmtPhotometry(w,bw,zp);
        end
        
    end
   
    enumeration
        % Photometric system
        %  Band ( wavelength , bandwidth , zero point )
        U  ( 0.360e-6 , 0.070e-6 , 2.0e12 )
        B  ( 0.440e-6 , 0.100e-6 , 5.4e12 )
        V  ( 0.550e-6 , 0.090e-6 , 3.3e12 )
        R  ( 0.640e-6 , 0.150e-6 , 4.0e12 )
        I  ( 0.790e-6 , 0.150e-6 , 2.7e12 )
        J  ( 1.215e-6 , 0.260e-6 , 1.9e12 )
        H  ( 1.654e-6 , 0.290e-6 , 1.1e12 )
        Ks ( 2.157e-6 , 0.320e-6 , 5.5e11 )
        K  ( 2.179e-6 , 0.410e-6 , 7.0e11 )
        L  ( 3.547e-6 , 0.570e-6 , 2.5e11 )
        M  ( 4.769e-6 , 0.450e-6 , 8.4e10 )
        Na ( 0.589e-6 , 0        , 3.3e12 )
    end
    
end