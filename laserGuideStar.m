classdef laserGuideStar < source
    %LASERGUIDESTAR Create a Laser Guide Star object
    %
    % lgs = laserGuideStar(apertureSize,apertureDistance,...
    %        meanAltitude,fwhmInArcsec,sourceParams) creates a Laser Guide
    %        Star object with the size of the aperture used to image the
    %        LGS, the distance between the LGS launch telescope and the
    %        furthest aperture, the LGS mean altitude, the LGS
    %        intensity profil FHWM in arcsec and the usual source
    %        parameters
    %
    % Example: A LGS on a 25m telescope with a 60x60 Shack-Hartmann WFS
    % launched from the telscope edge at an altitude of 90km and a fwhm of
    % 1 arcsec is created with:
    % lgs = laserGuideStar(25/60,25,90e3,1,...
    %        'wavelength',photometry.Na,...
    %        'height',1e3*((-5:5)+90),...
    %        'viewPoint',[-25/2,0]
    % The telescope pupil needs to be sampled with at leat 60*44 pixels
    %
    % See also: source
    
methods
    
    function obj = laserGuideStar(apertureSize,apertureDistance,...
            meanAltitude,fwhmInArcsec,varargin)
        
        obj = obj@source(varargin{:});
        
        if nargin>1
            set(obj,'objectiveFocalLength',meanAltitude)
            
            fwhm =  skyAngle( fwhmInArcsec, 'ARCSEC');
            width = skyAngle( apertureDistance*...
                (obj(1,1,end).height-obj(1,1,1).height)/meanAltitude^2 );
            fprintf(' LGS angular size: %5.3farcsec\n',width.ARCSEC);
            
            apertureImageFwhm = ...
                skyAngle( obj(1).wavelength/apertureSize );
            fprintf(' LGS aperture image FWHM: %5.3farcsec\n',apertureImageFwhm.ARCSEC)
            
            pixelWidth = ceil( 2*width/apertureImageFwhm );
            fprintf(' LGS size in pixel: %d\n',pixelWidth);
            
            if ~isempty(fwhm)
                pixelFwhm = 2*fwhm/apertureImageFwhm;
                fprintf(' LGS FWHM in pixel: %5.3f\n',pixelFwhm);
                intensityProfile = tools.gaussian(pixelWidth,pixelFwhm,...
                    ceil(pixelFwhm/(2*sqrt(2*log(2)))*3));
                set(obj,'extent',intensityProfile)
            end
            
        end
        
    end
    
end
    
end

