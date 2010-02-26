classdef influenceFunction < handle
    
    properties
        % mechanicalCoupling
        mechCoupling;
        % spline polynomials coefs
        splineP;
        % path listener
        bezierListener;
        % modes
        modes
        % influence function tag
        tag = 'BEZIER INFLUENCE FUN';
    end
    
    properties (SetObservable=true)
        % path
        bezier;
    end
    
    properties (Dependent)
        % points
        P;
    end
    
    properties (Access=private)
        % points
        p_P;
        % influence function display handle;
        displayHandle
        log
    end
    
    methods
        
        %% Constructor
        function obj = influenceFunction(points,mech)
            obj.p_P = zeros(7,2);
            obj.p_P(1,:) = [0,1];
            obj.p_P(2,2) = 1;
            obj.p_P(6,2) = 0;
            obj.p_P(7,:) = [2,0];
            if ischar(points)
                switch points
                    case 'overshoot'
                        points = {0.2,[0.4,0.7],[0.5,0.4],0.3,1};
                    case 'monotonic'
                        points = {0.2,[0.4,0.7],[0.6,0.4],1,1};
                end
            end
            obj.mechCoupling = mech;
            obj.P = points;
            obj.bezierListener = addlistener(obj,'bezier','PostSet',...
                @(src,evnt) obj.show );
            obj.bezierListener.Enabled = false;
            obj.log = logBook.checkIn(obj);
        end
        
        %% Destructor
        function delete(obj)
            if ishandle(obj.displayHandle)
                delete(get(obj.displayHandle,'parent'));
            end
            checkOut(obj.log,obj)
        end
        
        %% Set and Get P properties
        function out = get.P(obj)
            out = obj.p_P;
        end
        function set.P(obj,val)
            obj.p_P(2,1) = val{1};
            obj.p_P(3,:) = val{2};
            obj.p_P(4,:) = val{3};
            obj.p_P(5,:) = (-1/val{4})*obj.p_P(3,:)+(1+1/val{4})*obj.p_P(4,:);
            obj.p_P(6,1) = val{5};
            t = linspace(0,1,101)';
            obj.bezier = ...
                ((1-t).^3)*obj.p_P(1,:) + ...
                3.*((1-t).^2.*t)*obj.p_P(2,:) + ...
                3.*((1-t).*t.^2)*obj.p_P(3,:) + ...
                (t.^3)*obj.p_P(4,:);
            t(1) = [];
            obj.bezier = [obj.bezier ; ...
                ((1-t).^3)*obj.p_P(4,:) + ...
                3.*((1-t).^2.*t)*obj.p_P(5,:) + ...
                3.*((1-t).*t.^2)*obj.p_P(6,:) + ...
                (t.^3)*obj.p_P(7,:) ];
            obj.bezier(:,1) = obj.bezier(:,1)/...
                spline(obj.bezier(:,2),obj.bezier(:,1),obj.mechCoupling);
            u = [-flipud(obj.bezier(:,1)) ; obj.bezier(2:end,1)];
            v = [flipud(obj.bezier(:,2)) ; obj.bezier(2:end,2)];
            obj.splineP = spline(u,v);
            
        end
        
        function out = mtimes(obj,c)
            out = obj.modes*c;
        end
        function out = mldivide(obj,c)
            out = obj.modes\c;
        end
        
        function show(obj)
            if ishandle(obj.displayHandle)
                x = linspace(0,obj.bezier(end,1),101);
                set(obj.displayHandle(1),...
                    'XData',obj.bezier(:,1),...
                    'YData',obj.bezier(:,2))
                set(obj.displayHandle(2),...
                    'YData',ppval(obj.splineP,x))
            else
                x = linspace(0,obj.bezier(end,1),101);
                obj.displayHandle = plot(obj.bezier(:,1),obj.bezier(:,2),x,ppval(obj.splineP,x),'r--');
                grid
                xlabel('Normalized actuator pitch')
                ylabel('Normalized stroke')
            end
        end
        
        function setInfluenceFunction(obj,nIF,resolution,validActuator,ratioTelDm)
            if nargin<5
                ratioTelDm = 1;
            end
            z = linspace(-1,1,nIF)*(nIF-1)/2;
            u0 = ratioTelDm.*linspace(-1,1,resolution)*(nIF-1)/2;
            w = zeros(resolution,nIF);
            for kIF = 1:nIF
                u = u0 - z(kIF);
                index = u >= -obj.bezier(end,1) & u <= obj.bezier(end,1);
                w(index,kIF) = ppval(obj.splineP,u(index));
            end
            obj.modes = zeros(resolution^2,sum(validActuator(:)));
            kIF = 0;
            for jIF = 1:nIF
                for iIF = 1:nIF
                    if validActuator(iIF,jIF)
                        buffer = w(:,iIF)*w(:,jIF)';
                        kIF = kIF + 1;
                        obj.modes(:,kIF) = buffer(:);
                    end
                end
            end
        end
    end
    
end