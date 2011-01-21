classdef tomography < handle
    %% TOMOGRAPHY Create a tomography object
    %
    % tomo = tomography(sampling,diameter,atmModel,guideStar)
    % tomo = tomography(sampling,diameter,atmModel,guideStar,tomoStar)
    % tomo = tomography(...,'pupil',pupilMask,'unit',unit)
    
    properties
        tag = 'TOMOGRAPHY';
        % pupil sampling
        sampling
        % pupil diameter
        diameter
        % pupil mask
        pupil
        % atmosphere mode
        atmModel
        % guide stars
        guideStar
        % guide star covariance matrix
        Cxx
        % tomo star covariance matrix
        Coo
        % tomo star / guide star covariance matrix
        Cox
        nGuideStar
        nTomoStar
        % error unit
        unit;
        % Bayesian minimum mean square error        
        Bmse;
        % tip-tilt removal flag
        rmTipTilt
        % diameter of the LGS lauch telescope
        tipTiltDiameter
        % zernike / guide star covariance matrix
        Cax
        % zernike / tomo star covariance matrix
        Cao
        % zernike polynomials matrix
        Z
        % zernike covariance matrix
        Caa
        % Tomographic reconstructor
        tomoBuilder
    end
    
    properties (Dependent)
        % tomographic estimates directions
        tomoStar
    end
    
    properties (Dependent,SetAccess=private)
        % error variance
        var
        % error rms
        rms
        % error variance map
        varMap
        % error rms map
        rmsMap        
    end
    
    properties (Access=private)
        log
        p_tomoStar
    end
 
    methods
        
        %% Constructor
        function obj = tomography(sampling,diameter,atmModel,guideStar,varargin)
            
            inputs = inputParser;
            inputs.addRequired('sampling',@isnumeric);
            inputs.addRequired('diameter',@isnumeric);
            inputs.addRequired('atmModel',@(x) isa(x,'atmosphere'));
            inputs.addRequired('guideStar',@(x) isa(x,'source'));
            inputs.addOptional('tomoStar',[],@(x) isa(x,'source'));
            inputs.addParamValue('pupil',true(sampling),@islogical);
            inputs.addParamValue('unit',[],@isnumeric);
            inputs.addParamValue('rmTipTilt',false,@islogical);
            inputs.addParamValue('tipTiltDiameter',diameter,@isnumeric);
            
            inputs.parse(sampling,diameter,atmModel,guideStar,varargin{:});
            
            obj.sampling   = inputs.Results.sampling;
            obj.diameter   = inputs.Results.diameter;
            obj.atmModel   = inputs.Results.atmModel;
            obj.guideStar  = inputs.Results.guideStar;
            obj.nGuideStar = length(obj.guideStar);
            obj.p_tomoStar = inputs.Results.tomoStar;    
            obj.nTomoStar  = length(obj.p_tomoStar);
            obj.pupil      = inputs.Results.pupil;   
            obj.unit       = inputs.Results.unit;   
            obj.rmTipTilt  = inputs.Results.rmTipTilt;   
            obj.tipTiltDiameter  = inputs.Results.tipTiltDiameter;   
            
            obj.log = logBook.checkIn(obj);
            
            add(obj.log,obj,'Computing the covariance matrices')
            
            poolWasAlreadyOpen = true;
            if matlabpool('size')==0
                matlabpool('open')
                poolWasAlreadyOpen = false;
            end
            
            [obj.Cxx,obj.Cox] = ...
                phaseStats.spatioAngularCovarianceMatrix(...
                obj.sampling,obj.diameter,...
                obj.atmModel,obj.guideStar,obj.p_tomoStar,'mask',obj.pupil);
            obj.Coo = phaseStats.spatioAngularCovarianceMatrix(...
                obj.sampling,obj.diameter,...
                obj.atmModel,obj.p_tomoStar(1),'mask',obj.pupil);
            if obj.rmTipTilt
                obj.Cax = phaseStats.spatioAngularCovarianceMatrix(...
                    obj.sampling,obj.tipTiltDiameter,...
                    obj.atmModel,obj.guideStar,obj.guideStar,...
                    'mask',obj.pupil,'tipTilt',true);
                obj.Cao = phaseStats.spatioAngularCovarianceMatrix(...
                    obj.sampling,obj.tipTiltDiameter,...
                    obj.atmModel,obj.p_tomoStar,obj.guideStar,...
                    'mask',obj.pupil,'tipTilt',true);
                zern    = zernike(2:3,obj.tipTiltDiameter,'resolution',obj.sampling,'pupil',obj.pupil);
                obj.Z   = zern.p(obj.pupil,:);
                obj.Caa = phaseStats.zernikeAngularCovariance(zern,obj.atmModel,obj.guideStar);
            end
            
            add(obj.log,obj,'Computing the tomographic builder')
            
            m_tomoBuilder = cell(obj.nTomoStar,1);
            m_Cox = obj.Cox;
            m_Cxx = obj.Cxx;
            parfor k=1:obj.nTomoStar
                m_tomoBuilder{k} = m_Cox{k}/m_Cxx;
            end
            obj.tomoBuilder = m_tomoBuilder;
            
            if ~poolWasAlreadyOpen 
                matlabpool('close')
            end
            
        end
        
        %% Set/Get tomoStar
        function set.tomoStar(obj,val)
            obj.p_tomoStar = val;
            add(obj.log,obj,'Computing the tomo/guide stars covariance matrix')
            obj.Cox = phaseStats.spatioAngularCovarianceMatrix(...
                obj.sampling,obj.diameter,...
                obj.atmModel,obj.guideStar,obj.p_tomoStar,'mask',obj.pupil);
            obj.Coo = phaseStats.spatioAngularCovarianceMatrix(...
                obj.sampling,obj.diameter,...
                obj.atmModel,obj.p_tomoStar(1),'mask',obj.pupil);
        end
        function val = get.tomoStar(obj)
            val = obj.p_tomoStar;
        end
        
        function lmmse(obj,fieldStar)
            %% LMMSE Linear Minimum Mean Square Error
            
            add(obj.log,obj,'Bayesian Minimum Mean Square Error computation in progress...')
            
            if obj.rmTipTilt
                m_Cox = obj.Cox{1};
                Czx = cell2mat( cellfun(@(x) obj.Z*x , obj.Cax , 'uniformOutput' , false) );
                Czo = cell2mat( cellfun(@(x) obj.Z*x , obj.Cao, 'uniformOutput' , false) );
                Czz = cell2mat( cellfun(@(x) obj.Z*x*obj.Z',obj.Caa, 'uniformOutput' , false) );
                M = m_Cox/obj.Cxx;
                m_Bmse  = obj.Coo - M*m_Cox';
                term1 = M*Czo;
                term2 = M*(Czz - Czx - Czx')*M';
                obj.Bmse = { m_Bmse + term1 + term1' + term2 };
            else
                if nargin==1
                    fun = @(x,y) obj.Coo - y*x';
                    obj.Bmse = cellfun( fun , obj.Cox , obj.tomoBuilder , 'uniformOutput' , false );
                else
                    Co1x = ...
                        phaseStats.spatioAngularCovarianceMatrix(...
                        obj.sampling,obj.diameter,...
                        obj.atmModel,obj.guideStar,fieldStar,'mask',obj.pupil);
                    fun = @(x) obj.Coo + ...
                        obj.Cox{1}/obj.Cxx*obj.Cox{1}' - ...
                        x/obj.Cxx*obj.Cox{1}' - obj.Cox{1}/obj.Cxx*x';
                    obj.Bmse = cellfun( fun , Co1x , 'uniformOutput' , false );
                end
            end
            
        end
       
        function map = get.varMap(obj)
            %% VARMAP Pupil error variance map
            
            if isempty(obj.Bmse)
                lmmse(obj)
            end
            add(obj.log,obj,'Pupil error variance map computation in progress...')
            out = cellfun( @diag , obj.Bmse' , 'uniformOutput', false );
            out = cell2mat(out);
            map = zeros(obj.sampling,obj.sampling*obj.nTomoStar);
            mask = repmat( obj.pupil, 1 , obj.nTomoStar);
            map(mask) = out;
                
        end
       
        function out = get.var(obj)
            %% VAR Pupil error variance
            
            if isempty(obj.Bmse)
                lmmse(obj)
            end
            add(obj.log,obj,'Pupil error variance computation in progress...')
            if iscell(obj.Bmse)
                out = cellfun( @trace , obj.Bmse , 'uniformOutput', true );
            else
                out = trace(obj.Bmse);
            end
            out = out/sum(obj.pupil(:));
        end
        
        function map = get.rmsMap(obj)
            %% RMSMAP Pupil error rms map
            
            map = opd(obj,obj.varMap);
        end
        
        function out = get.rms(obj)
            %% RMS Pupil error rms
            
            out = opd(obj,obj.var);
        end
     
    end
       
    
    methods (Access=private)
        
        function out = opd(obj,val)
            out = sqrt(val);
            if ~isempty(obj.unit)
                add(obj.log,obj,sprintf('Rms in 1E%d meter',obj.unit))
                out = 10^-obj.unit*...
                    out*obj.atmModel.wavelength/(2*pi);
            end
        end
        
    end
    
end