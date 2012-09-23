classdef photometry
    properties
        wavelength % [micron]
        bandwidth  % [micron]
        zeroPoint  % [ph/s]
    end
    
    methods
        
        %% Constructor
        function obj = photometry(w,bw,zp)
            obj.wavelength = w;
            obj.bandwidth  = bw;
            obj.zeroPoint  = zp/368;
        end
        
        %% Get nPhoton
        function out = nPhoton(obj,magnitude)
            out = obj.zeroPoint*10^(-0.4*magnitude);
        end
        %% Get magnitude
        function out = magnitude(obj,nPhoton)
            out = -2.5*log10(nPhoton/obj.zeroPoint);
        end
        
        function out = rdivide(obj,val)
            %% WAVELENGTH DIVISION
            %
            % out = obj./val divides the wavelength by val
            
            out = obj.wavelength./val;
        end
        
        function out = mrdivide(obj,val)
            %% WAVELENGTH DIVISION
            %
            % out = obj/val divides the wavelength by val
            
            out = rdivide(obj,val);
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
        EOS ( 1.064e-6 , 0        , 3.3e12 )
        HK  ( (1.654e-6*1.1e12+2.179e-6*7.0e11)/(1.1e12+7.0e11) , 0.290e-6+00.410e-6 , 1.1e12+7.0e11 )
    end
    
end