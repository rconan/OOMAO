classdef (Sealed) logBook < handle
    
    properties 
        verbose=true;
        nLogs = 1000;
        queue;
        queueIndex;
        nObj = 0; 
        writable   = true;
    end
   
    properties (SetAccess=private)
        logs;
    end
    
    properties (Access=private)
        kLogs = 0;
    end
    
    methods (Access=private)

        % Constructor
        function obj = logBook(varargin)
            for k=1:nargin
                varargin{k}.logsId = k;
            end
            obj.queue = timer;
            obj.queue.BusyMode = 'queue';
%             obj.queue.StartDelay = 1;
            obj.queue.timerFcn = @(src,event) cellfun(@(x) fprintf(x), obj.logs(obj.queueIndex));
            obj.queue.stopFcn = @stopQueue;
            function stopQueue(src,event)
                obj.queueIndex = [];
            end
        end
        
    end
    
    methods
        
        function delete(obj)
            delete(obj.queue)
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
        
        % Objects info. listener
        function add(obj,src,info)
            k           = mod(obj.kLogs,obj.nLogs) + 1;
%             obj.logs{k} = sprintf(' @(%s)> %s\n',src.logsId,evnt.info);
            obj.logs{k} = sprintf(' @(%s)> %s\n',lower(src.tag),info);
            if obj.verbose
%                 fprintf(obj.logs{k});
                obj.queueIndex = [obj.queueIndex, k];
                if strcmp(obj.queue.Running,'off')
                    start(obj.queue);
                end
            end
            obj.kLogs = obj.kLogs + 1;
        end 
        
        function checkOut(obj,src)
            if obj.verbose
                add(obj,src,'Terminated!')
            end
            obj.nObj = obj.nObj - 1;
            if obj.nObj==0
                wait(obj.queue)
                delete(obj)
                fprintf(' @(logBook)> Closing the log book!\n')
                fprintf('~~~~~~~~~~~~~~~~~~~\n OOMAO''S GONE!\n~~~~~~~~~~~~~~~~~~~\n')
            end
        end

    end
    
    methods (Static)
        
        function obj = checkIn(src)
            persistent localObj
            if isempty(localObj) || ~isvalid(localObj)
                fprintf('~~~~~~~~~~~~~~~~~~~\n BEWARE OF OOMAO!\n~~~~~~~~~~~~~~~~~~~\n')
                fprintf(' @(logBook)> Opening the log book!\n')
                localObj = logBook;
            end
            obj = localObj;
            if nargin>0 && obj.writable%&& isvalid(src)
                nSrc = numel(src);
                obj.nObj = obj.nObj + nSrc;
                if nSrc>1
                    add(obj,src(1),sprintf('%d created!',nSrc))
                else
                    add(obj,src,'Created!')
                end
            end
        end
        
        function PAUSE
            log = logBook.checkIn;
            log.writable = false;
        end
        
        function RESUME
            log = logBook.checkIn;
            log.writable = true;
        end
      
%         function add(fromObj,info)
% %             if all(fromObj.logsId<0) % Check if object is already referenced in the logBook
% %                 % Look for the object in the base workspace
% %                 s = evalin('base','whos');
% %                 % Check if the logBook exists in the base workspace
% %                 index = strmatch('logBook',{s.class});
% %                 if isempty(index) || evalin('base','~isvalid(cougarLogs)')
% %                     assignin('base','cougarLogs',logBook);
% %                     fprintf('~~~~~~~~~~~~~~~~~~~\n BEWARE OF COUGAR!\n~~~~~~~~~~~~~~~~~~~\n')
% %                 end
% %                 % Check if fromObj class exists in the base workspace
% %                 index = strmatch(class(fromObj),{s.class})';
% %                 % if not do nothing and go back to the calling function
% %                 if isempty(index)
% %                     return
% %                 else
% %                     % if the object class is found, go for the class
% %                     % object
% %                     fromObj.logsId = Inf;
% %                     for k=index
% %                         flag =  evalin('base',...
% %                             ['isvalid(',s(k).name,') && all(isinf(',s(k).name,'.logsId))']);
% %                         if flag
% %                             break;
% %                         end
% %                     end
% %                     % if the class object is not found, go back to the calling function
% %                     if ~flag
% %                         fromObj.logsId = -1;
% %                         return
% %                     end
% %                     % if it is found, set its Id and sign in it to the logbook
% %                     fromObj.logsId = s(k).name;
% %                     evalin('base',['cougarLogs.signIn(',s(k).name,');']); % Sign in the object in the logBook
% %                 end
% %             end
% %             t = timerfind('name','logSignIn');
% %             if ~isempty(t)
% %                 wait(t)
% %                 delete(t)
% %             end
%             % Queue the object message
%             notify(fromObj,'logs',logBookData(info));
%         end
% 
%         function signIn(fromObj)
% 
% %             start(timer('BusyMode','queue','StartDelay',1,'name','logSignIn',...
% %                 'TimerFcn',{@timerCallBack, fromObj}))%,...
% %            %     'StopFcn','delete(timerfind(''name'',''logSignIn''))'))
% %             function timerCallBack( timerObj, event, fromObj)
% % %                                 logBook.add(fromObj,'Started!' )
% %             if all(fromObj.logsId<0) % Check if object is already referenced in the logBook
%                 % Look for cougar in the base workspace
%                 s = evalin('base','whos');
%                 % Check if the logBook exists in the base workspace
%                 index = strmatch('logBook',{s.class});
%                 if isempty(index) || evalin('base','~isvalid(cougarLogs)')
%                     evalin('base','global cougarLogs');
%                     assignin('base','cougarLogs',logBook);
%                     fprintf('~~~~~~~~~~~~~~~~~~~\n BEWARE OF COUGAR!\n~~~~~~~~~~~~~~~~~~~\n')
%                 end
%                 global cougarLogs
%                 cougarLogs.newEntry(fromObj);
% % %                 % Look for the object in the base workspace
% % %                 s = evalin('caller','whos');
% %                 % Check if fromObj class exists in the base workspace
% %                 index = strmatch(class(fromObj),{s.class})';
% %                 % if not do nothing and go back to the calling function
% %                 if isempty(index)
% %                     return
% %                 else
% %                     % if the object class is found, go for the class
% %                     % object
% %                     fromObj.logsId = Inf;
% %                     for k=index
% %                         flag =  evalin('base',...
% %                             ['isvalid(',s(k).name,') && all(isinf(',s(k).name,'.logsId))']);
% %                         if flag
% %                             break;
% %                         end
% %                     end
% %                     % if the class object is not found, go back to the calling function
% %                     if ~flag
% %                         fromObj.logsId = -1;
% %                         return
% %                     end
% %                     % if it is found, set its Id and sign in it to the logbook
% %                     fromObj.logsId = s(k).name;
% %                     evalin('base',['cougarLogs.newEntry(',s(k).name,');']); % Sign in the object in the logBook
% %                 end
% %             end
%             % Queue the object message
%             notify(fromObj,'logs',logBookData('Started!'));           
% %             end
%         end
    end

end

