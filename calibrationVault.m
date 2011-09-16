classdef calibrationVault < handle
    %% CALIBRATIONVAULT Create a calibrationVault object
    %
    % calib = calibrationVault(dm,wfs,calibMatrix) create a
    % calibrationVault object storing the calibration matrix between the dm
    % and the wfs. The svd of calibMatrix is computed and used to compute
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
        eigLine;
        p_threshold;
        p_nThresholded;
    end
    
    methods
        
        %% Constructor
        function obj = calibrationVault(calibMatrix)
            
            obj.D      = calibMatrix;
            obj.log    = logBook.checkIn(obj);
            
            add(obj.log,obj,'Computing the SVD of the calibration matrix!')
            
            [obj.U,S,obj.V] = svd(calibMatrix,0);
            obj.eigenValues = diag(S);
            
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
            subplot(1,2,1)
            imagesc(obj.D)
            xlabel('DM actuators')
            ylabel('WFS slopes')
            ylabel(colorbar,'slopes/actuator stroke')
            obj.eigAxis = subplot(1,2,2);
            semilogy(obj.eigenValues,'.')
            xlabel('Eigen modes')
            ylabel('Eigen values')
            
            if ~isempty(obj.eigLine)
                obj.eigLine = line(get(obj.eigAxis,'xlim'),ones(1,2)*obj.p_threshold,'color','r','parent',obj.eigAxis);
            end
            drawnow

        end
        
    end
    
    methods (Access=private)
        
        function updateCommandMatrix(obj)
            %% UPDATECOMMANDMATRIX Update the command matrix
            
            figure(get(obj.eigAxis,'parent'))
            if isempty(obj.eigLine)
                obj.eigLine = line(get(obj.eigAxis,'xlim'),ones(1,2)*obj.p_threshold,'color','r','parent',obj.eigAxis);
            else
                set(obj.eigLine,'ydata',ones(1,2)*obj.p_threshold)
            end
            drawnow

            add(obj.log,obj,'Updating the command matrix!')
    
            nEigenValues = length(obj.eigenValues) - obj.nThresholded;
            u = 1:nEigenValues;
            iS = diag(1./obj.eigenValues(u));
            obj.M = obj.V(:,u)*iS*obj.U(:,u)';
            obj.M = obj.spaceJump*obj.M;
        end
        
    end
    
    methods (Static)
        
        function obj = loadobj(obj)
            show(obj)
        end
    end
end