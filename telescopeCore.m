classdef telescopeCore < handle
    % TELESCOPECORE Create a telescopeCore object
    %
    % tel = telescopeCore(D) creates a telescopeCore object from the diameter D.
    %
    % tel = telescopeCore(D,'parameter',value,...) creates a telescopeCore object from
    % the diameter D and from optionnal parameter-value pair arguments. The
    % parameters are obstructionRatio, fieldOfViewInArcsec, fieldOfViewInArcmin
    % or resolution.
    
    properties
        % diameter
        D;
        % central obstruction ratio
        obstructionRatio;
        % conjugation altitude
        conjugationHeight;
        % focalisation distance
        focalDistance;
        % field-of-view
        fieldOfView;
        % diameter resolution in pixel
        resolution;
    end
    
    properties (Dependent)
        % telescope pupil mask
        pupil;
    end
    
    properties (Dependent,SetAccess=private)
        % radius
        R;
        % telescope pupil mask
        pupilLogical;
        % telescope area
        area;
    end
    
    properties (Access=protected)
        p_pupil;
    end
    
    methods
        
        % Constructor
        function obj = telescopeCore(D,varargin)
            p = inputParser;
            p.addRequired('D', @isnumeric);
            p.addParamValue('obstructionRatio', 0, @isnumeric);
            p.addParamValue('conjugationHeight', 0, @isnumeric);
            p.addParamValue('focalDistance', Inf, @isnumeric);
            p.addParamValue('fieldOfViewInArcsec', [], @isnumeric);
            p.addParamValue('fieldOfViewInArcmin', [], @isnumeric);
            p.addParamValue('resolution', [], @isnumeric);
            p.parse(D, varargin{:});
            obj.D                = p.Results.D;
            obj.conjugationHeight = p.Results.conjugationHeight;
            obj.focalDistance = p.Results.focalDistance;
            obj.obstructionRatio = p.Results.obstructionRatio;
            if ~isempty(p.Results.fieldOfViewInArcsec)
                obj.fieldOfView      = p.Results.fieldOfViewInArcsec./cougarConstants.radian2arcsec;
            elseif ~isempty(p.Results.fieldOfViewInArcmin)
                obj.fieldOfView      = p.Results.fieldOfViewInArcmin./cougarConstants.radian2arcmin;
            else
                obj.fieldOfView      = 0;
            end
            obj.resolution       = p.Results.resolution;
        end

        % Get and Set the pupil
        function pupil = get.pupil(obj)
            pupil = obj.p_pupil;
            if isempty(pupil) && ~isempty(obj.resolution)
                pupil = utilities.piston(obj.resolution);
                if obj.obstructionRatio>0
                    pupil = pupil - ...
                        utilities.piston(...
                        round(obj.resolution.*obj.obstructionRatio),...
                        obj.resolution);
                end
            end
        end
        function set.pupil(obj,val)
            obj.p_pupil = val;
        end
        
        function pupilLogical = get.pupilLogical(obj)
            pupilLogical = logical(obj.pupil>0);
        end
        
        function out = get.R(obj)
            out = obj.D/2;
        end
        
        function out = get.area(obj)
            out = pi*obj.R^2*(1-obj.obstructionRatio^2);
        end
        
        function out = diameterAt(obj,height)
            out = obj.D + 2.*height.*tan(obj.fieldOfView/2);
        end
        
        function out = getWave(obj,src)
            % GETWAVE Wave propagation to the telescope
            %
            % out = getWave(obj,src) computes the propagation from the
            % source to the telescope pupil
            
            srcDims = size(src);
            
            out = bsxfun( @times , repmat( obj.pupil ,[1 srcDims(2)/srcDims(1)] ) ,...
                exp( 1i.*cell2mat( ...
                cellfun(@(x) fresnelPropagation(x,obj),num2cell(src),...
                'UniformOutput',false) ) ) );
        end
        
        function out = otf(obj, r)
            % OTF Telescope optical transfert function
            %
            % out = otf(obj, r) Computes the telescope optical transfert function
            
%             out = zeros(size(r));
            if obj.obstructionRatio ~= 0
                out = pupAutoCorr(obj.D) + pupAutoCorr(obj.obstructionRatio*obj.D) - ...
                    2.*pupCrossCorr(obj.D./2,obj.obstructionRatio*obj.D./2);
            else
                out = pupAutoCorr(obj.D);
            end
            out = out./(pi*obj.D*obj.D.*(1-obj.obstructionRatio*obj.obstructionRatio)./4);
            
            function out1 = pupAutoCorr(D)
                
                index       = r <= D;
                red         = r(index)./D;
                out1        = zeros(size(r));
                out1(index) = D.*D.*(acos(red)-red.*sqrt((1-red.*red)))./2;
                
            end
            
            
            function out2 = pupCrossCorr(R1,R2)
                
                out2 = zeros(size(r));
                
                index       = r <= abs(R1-R2);
                out2(index) = pi*min([R1,R2]).^2;
                
                index       = (r > abs(R1-R2)) & (r < (R1+R2));
                rho         = r(index);
                red         = (R1*R1-R2*R2+rho.*rho)./(2.*rho)/(R1);
                out2(index) = out2(index) + R1.*R1.*(acos(red)-red.*sqrt((1-red.*red)));
                red         = (R2*R2-R1*R1+rho.*rho)./(2.*rho)/(R2);
                out2(index) = out2(index) + R2.*R2.*(acos(red)-red.*sqrt((1-red.*red)));
                
            end
            
        end
        
        function out = psf(obj,f)
            % PSF Telescope point spread function
            %
            % out = psf(obj, f) computes the telescope point spread function
            
            out   = ones(size(f)).*pi.*obj.D.^2.*(1-obj.obstructionRatio.^2)./4;
            index = find(f);
            u = pi.*obj.D.*f(index);
            surface = pi.*obj.D.^2./4;
            out(index) = surface.*2.*besselj(1,u)./u;
            if obj.obstructionRatio>0
                u = pi.*obj.D.*obj.obstructionRatio.*f(index);
                surface = surface.*obj.obstructionRatio.^2;
                out(index) = out(index) - surface.*2.*besselj(1,u)./u;
                
            end
            out = abs(out).^2./(pi.*obj.D.^2.*(1-obj.obstructionRatio.^2)./4);
        end
        
        
        function out = fullWidthHalfMax(obj)
            % FULLWIDTHHALFMAX Full Width at Half the Maximum evaluation
            %
            % out = fullWidthHalfMax(a) computes the FWHM of a telescope
            % object. Units are m^{-1}. To convert it in arcsecond,
            % multiply by the wavelength then by radian2arcsec.
            
            x0 = [0,2/obj.D];
            [x,fval,exitflag] = fzero(@(x) psf(obj,abs(x)./2) - psf(obj,0)./2,x0,optimset('TolX',1e-9));
            if exitflag<0
                warning('cougar:telescope:fullWidthHalfMax',...
                    'No interval was found with a sign change, or a NaN or Inf function value was encountered during search for an interval containing a sign change, or a complex function value was encountered during the search for an interval containing a sign change.')
            end
            out = abs(x);
        end
        
        function out = entrappedEnergy(obj,eHalfSize,trap,psfOrOtf)
            % ENTRAPPEDENERGY Encircled of ensquared energy
            %
            % out = entrappedEnergy(obj,eHalfSize,trap) computes the
            % entraped energy in a circle of radius eHalfSize if trap is
            % set to 'circle' or in a square of half length eHalfSize if
            % trap is set to 'square'
            
            if nargin<4
                psfOrOtf = 'psf';
            end
            switch lower(psfOrOtf)
                case 'otf'
                     switch lower(trap)
                        case 'circle'
                            out = 2*pi*quadgk( ...
                                @(r) r.*otf(obj,r).*...
                                (2.*besselj(1,2*pi.*eHalfSize.*r)./(2*pi.*eHalfSize.*r)),0,obj.D)*...
                                pi*eHalfSize^2;
                        case 'square'
                            a = 2*eHalfSize;
                            out = quad2d(...
                                @(o,r) r.*otf(obj,r).*...
                                (sin(pi.*r.*cos(o).*a)./(pi.*r.*cos(o).*a)).*...
                                (sin(pi.*r.*sin(o).*a)./(pi.*r.*sin(o).*a)), ...
                                0,2*pi,0,obj.D).*a.*a;
                        otherwise
                            error('cougar:telescope:entrapedEnergy',...
                                'The trap is either a circle or a square!')
                    end
               otherwise
                    switch lower(trap)
                        case 'circle'
                            out = quadgk(@(x)x.*psf(obj,x),0,eHalfSize)*2*pi;
                        case 'square'
                            out = quad2d(@(x,y)psf(obj,hypot(x,y)),0,eHalfSize,0,eHalfSize)*4;
                        otherwise
                            error('cougar:telescope:entrapedEnergy',...
                                'The trap is either a circle or a square!')
                    end
            end
        end
                
    end
    
end
