classdef (Hidden=true) logBook < handle
    
    properties 
        verbose=true;
        nLogs = 1000;
        queue;
        queueIndex;
    end
   
    properties (SetAccess=private)
        logs;
   end
   
    properties (SetAccess=private,GetAccess=private)
        hlInfo
        kLogs = 0;
    end
    
    methods

        % Constructor
        function obj = logBook(varargin)
            for k=1:nargin
                varargin{k}.logsId = k;
            end
            obj.hlInfo = event.listener(varargin,'logs',...
                @(src,evnt)listenInfo(obj,src,evnt));
            obj.queue = timer;
            obj.queue.BusyMode = 'queue';
%             obj.queue.StartDelay = 1;
            obj.queue.timerFcn = @(src,event) cellfun(@(x) fprintf(x), obj.logs(obj.queueIndex));
            obj.queue.stopFcn = @stopQueue;
            function stopQueue(src,event)
                obj.queueIndex = [];
            end
        end
        
        function delete(obj)
            delete(obj.queue)
        end
        
        % Add a new object to listener
        function newEntry(obj,newObj)
            if isobject(obj.hlInfo) % Check if it is already a listener
                % Looking for invalid objects
                invalidObj = cellfun(@(x) ~isvalid(x),obj.hlInfo.Source);
                if any(invalidObj)
                    index     = find(invalidObj);
                    nSignedIn = index(1);
                else
                    nSignedIn = length(obj.hlInfo.Source) + 1;
                end
            else
                nSignedIn = 1;
            end
%             newObj.logsId = nSignedIn;
            obj.hlInfo.Source{nSignedIn} = newObj;
            obj.hlInfo = event.listener(obj.hlInfo.Source,'logs',...
                @(src,evnt)listenInfo(obj,src,evnt));
        end
        
        % Display the last logs
        function tail(obj,n)
            if nargin<2
                a = 1;
            else
                a = max(length(obj.logs)-n+1,1);
            end
            cellfun(@(x) fprintf(x), obj.logs(a:end) );
        end
        
        % Resume logging
        function resume(obj)
             obj.hlInfo.Enabled = true;
        end
        
        % Pause logging
        function pause(obj)
             obj.hlInfo.Enabled = false;
        end

        % Objects info. listener
        function listenInfo(obj,src,evnt)
            k           = mod(obj.kLogs,obj.nLogs) + 1;
%             obj.logs{k} = sprintf(' @(%s)> %s\n',src.logsId,evnt.info);
            obj.logs{k} = sprintf(' @(%s)> %s\n',class(src),evnt.info);
            if obj.verbose
%                 fprintf(obj.logs{k});
                obj.queueIndex = [obj.queueIndex, k];
                if strcmp(obj.queue.Running,'off')
                    start(obj.queue);
                end
            end
            obj.kLogs = obj.kLogs + 1;
            checkWhatsLeft(obj)
        end 
        
        function checkWhatsLeft(obj)
            invalidObj = cellfun(@(x) ~isvalid(x),obj.hlInfo.Source);
            if all(invalidObj)
                wait(obj.queue)
                delete(obj)
                fprintf('~~~~~~~~~~~~~~~~~~~\n COUGAR''S GONE!\n~~~~~~~~~~~~~~~~~~~\n')
            end
        end

    end
    
    methods (Static)
        function add(fromObj,info)
%             if all(fromObj.logsId<0) % Check if object is already referenced in the logBook
%                 % Look for the object in the base workspace
%                 s = evalin('base','whos');
%                 % Check if the logBook exists in the base workspace
%                 index = strmatch('logBook',{s.class});
%                 if isempty(index) || evalin('base','~isvalid(cougarLogs)')
%                     assignin('base','cougarLogs',logBook);
%                     fprintf('~~~~~~~~~~~~~~~~~~~\n BEWARE OF COUGAR!\n~~~~~~~~~~~~~~~~~~~\n')
%                 end
%                 % Check if fromObj class exists in the base workspace
%                 index = strmatch(class(fromObj),{s.class})';
%                 % if not do nothing and go back to the calling function
%                 if isempty(index)
%                     return
%                 else
%                     % if the object class is found, go for the class
%                     % object
%                     fromObj.logsId = Inf;
%                     for k=index
%                         flag =  evalin('base',...
%                             ['isvalid(',s(k).name,') && all(isinf(',s(k).name,'.logsId))']);
%                         if flag
%                             break;
%                         end
%                     end
%                     % if the class object is not found, go back to the calling function
%                     if ~flag
%                         fromObj.logsId = -1;
%                         return
%                     end
%                     % if it is found, set its Id and sign in it to the logbook
%                     fromObj.logsId = s(k).name;
%                     evalin('base',['cougarLogs.signIn(',s(k).name,');']); % Sign in the object in the logBook
%                 end
%             end
%             t = timerfind('name','logSignIn');
%             if ~isempty(t)
%                 wait(t)
%                 delete(t)
%             end
            % Queue the object message
            notify(fromObj,'logs',logBookData(info));
        end

        function signIn(fromObj)

%             start(timer('BusyMode','queue','StartDelay',1,'name','logSignIn',...
%                 'TimerFcn',{@timerCallBack, fromObj}))%,...
%            %     'StopFcn','delete(timerfind(''name'',''logSignIn''))'))
%             function timerCallBack( timerObj, event, fromObj)
% %                                 logBook.add(fromObj,'Started!' )
%             if all(fromObj.logsId<0) % Check if object is already referenced in the logBook
                % Look for cougar in the base workspace
                s = evalin('base','whos');
                % Check if the logBook exists in the base workspace
                index = strmatch('logBook',{s.class});
                if isempty(index) || evalin('base','~isvalid(cougarLogs)')
                    evalin('base','global cougarLogs');
                    assignin('base','cougarLogs',logBook);
                    fprintf('~~~~~~~~~~~~~~~~~~~\n BEWARE OF COUGAR!\n~~~~~~~~~~~~~~~~~~~\n')
                end
                global cougarLogs
                cougarLogs.newEntry(fromObj);
% %                 % Look for the object in the base workspace
% %                 s = evalin('caller','whos');
%                 % Check if fromObj class exists in the base workspace
%                 index = strmatch(class(fromObj),{s.class})';
%                 % if not do nothing and go back to the calling function
%                 if isempty(index)
%                     return
%                 else
%                     % if the object class is found, go for the class
%                     % object
%                     fromObj.logsId = Inf;
%                     for k=index
%                         flag =  evalin('base',...
%                             ['isvalid(',s(k).name,') && all(isinf(',s(k).name,'.logsId))']);
%                         if flag
%                             break;
%                         end
%                     end
%                     % if the class object is not found, go back to the calling function
%                     if ~flag
%                         fromObj.logsId = -1;
%                         return
%                     end
%                     % if it is found, set its Id and sign in it to the logbook
%                     fromObj.logsId = s(k).name;
%                     evalin('base',['cougarLogs.newEntry(',s(k).name,');']); % Sign in the object in the logBook
%                 end
%             end
            % Queue the object message
            notify(fromObj,'logs',logBookData('Started!'));           
%             end
        end
    end

end

