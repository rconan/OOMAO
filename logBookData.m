classdef logBookData < event.EventData
   
    properties
        info;
    end
    
    methods

        function data = logBookData(newInfo)
            data.info = newInfo;
        end
        
    end
end