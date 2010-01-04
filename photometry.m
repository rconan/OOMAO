classdef photometry
% PHOTOMETRY Photometric system in astronomy

    properties (Constant=true)
        V =  550e-9;
        R =  700e-9;
        I =  900e-9;
        J = 1250e-9;
        H = 1650e-9
        K = 2200e-9
        L = 3400e-9
        M = 5000e-9;
    end
    properties (Constant=true,Hidden)
        % W/m^2/micron
        Ve0 = 3.92e-08;
        Re0 = 1.76e-08;
        Ie0 = 8.30e-09
        Je0 = 3.40e-09
        He0 = 7.00e-10
        Ke0 = 3.90e-10;
        Le0 = 8.10e-11
        Me0 = 2.20e-11;
        wavelengths = [550 700 900 1250 1650 2200 3400 5000]*1e-9;
        deltaWavelengths = [89 220 240 300 350 400 550 300]*1e-9;
        e0 = [3.92e-08 1.76e-08 8.30e-09 3.40e-09 7.00e-10 3.90e-10 8.10e-11 2.20e-11];
    end
end