classdef calibrationVault < handle
    %% CALIBRATIONVAULT Create a calibrationVault object
    %
    % calib = calibrationVault(calibMatrix) create a
    % calibrationVault object storing the calibration matrix between a dm
    % and a wfs. The svd of calibMatrix is computed and used to compute
    % the dm command matix
    
    properties
        % the calibration matrix
        D;
        % the SVD decomposition of the calibration matrix
        U;
        eigenValues;
        V;
        % the command matrix
        M;
        % the truncated calibration matrix based on the threshold of SVD eigen values
        truncD;
        % a matrix to project the command matrix in another sub-space
        spaceJump = 1;
        % tag
        tag = 'CALIBRATION VAULT';
        % DM influence function
        modes;
    end
    
    properties (Dependent)
        % the SVD threshold
        threshold;
        % the number of tresholded eigen values
        nThresholded;
    end
    
    properties (Access=private)
        log;
        eigAxis;
        eigModeAxis;
        eigLine;
        eigImage;
        p_threshold;
        p_nThresholded;
    end
    
    methods
        
        %% Constructor
        function obj = calibrationVault(calibMatrix,modes,pupil)
            
            obj.D      = calibMatrix;
            if nargin>1
                obj.modes  = bsxfun( @times, modes, pupil(:) );
            end
            obj.log    = logBook.checkIn(obj);
            
            add(obj.log,obj,'Computing the SVD of the calibration matrix!')
            
            [obj.U,S,obj.V] = svd(calibMatrix,0);
            obj.eigenValues = diag(S);
            
            iS = diag(1./obj.eigenValues);
            obj.M = obj.V*iS*obj.U';
            obj.p_nThresholded = 0;
            obj.p_threshold = obj.eigenValues(end);
            
            show(obj)
            
        end
        
        %% Destructor
        function delete(obj)
            checkOut(obj.log,obj)
        end
        
        
        %% Set/Get threshold
        function set.threshold(obj,val)
            obj.p_threshold = val;
            obj.p_nThresholded = sum(obj.eigenValues<val);
            updateCommandMatrix(obj)
        end
        function val = get.threshold(obj)
            val = obj.p_threshold;
        end
        
        %% Set/Get nTthresholded
        function set.nThresholded(obj,val)
            obj.p_nThresholded = val;
            obj.p_threshold = obj.eigenValues(end-val);
            updateCommandMatrix(obj)
        end
        function val = get.nThresholded(obj)
            val = obj.p_nThresholded;
        end
        
          
        function show(obj)
            
            figure
            
            subplot(2,2,[1,3])
            imagesc(obj.D)
            xlabel('DM actuators')
            ylabel('WFS slopes')
            ylabel(colorbar,'slopes/actuator stroke')
            
            obj.eigAxis = subplot(2,2,2);
            semilogy(obj.eigenValues,'.')
            xlabel('Eigen modes')
            ylabel('Eigen values')
            
            obj.eigLine = line(get(obj.eigAxis,'xlim'),ones(1,2)*obj.p_threshold,'color','r','parent',obj.eigAxis);
            
            if ~isempty(obj.modes)
                obj.eigModeAxis = subplot(2,2,4);
                obj.eigImage = imagesc(tools.toggleFrame(obj.modes*obj.V(:,end-obj.p_nThresholded)));
                axis square
                colorbar
            end
            
            drawnow

        end
        
    end
    
    methods (Access=private)
        
        function updateCommandMatrix(obj)
            %% UPDATECOMMANDMATRIX Update the command matrix
            
            figure(get(obj.eigAxis,'parent'))
%             if isempty(obj.eigLine)
%                 obj.eigLine = line(get(obj.eigAxis,'xlim'),ones(1,2)*obj.p_threshold,'color','r','parent',obj.eigAxis);
%             else
                set(obj.eigLine,'ydata',ones(1,2)*obj.p_threshold)
                if ~isempty(obj.modes)
                    set(obj.eigImage,'Cdata',tools.toggleFrame(obj.modes*obj.V(:,end-obj.p_nThresholded)))
                end
%             end
            drawnow

            add(obj.log,obj,'Updating the command matrix!')
    
            nEigenValues = length(obj.eigenValues) - obj.nThresholded;
            u = 1:nEigenValues;
            iS = diag(1./obj.eigenValues(u));
            obj.M = obj.V(:,u)*iS*obj.U(:,u)';
            obj.M = obj.spaceJump*obj.M;
            
            obj.truncD = obj.U(:,u)*diag(obj.eigenValues(u))*obj.V(:,u)';
        end
        
    end
    
    methods (Static)
        
        function obj = loadobj(obj)
            show(obj)
        end
    end
end