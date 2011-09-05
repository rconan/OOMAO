classdef sourceorama < handle
   
    properties
        
        % source object array
        srcs;
        % telescope+atmosphere object
        tel;
        % HDF5 file name
        hdf5file;
        % whole data time length
        sequenceTimeLength
        % tag
        tag = 'SOURCEORAMA';
    end
    
    properties (Dependent)
        % time length of each data set [s]
        dataSetTimeLength;
        
    end
    
    properties (SetAccess=private)
        % memory size of each data set [GB]
        dataSetMemSize;
    end
    
    properties (Access=private)
        % time length of each data set [s]
        p_dataSetTimeLength;
        log;
    end
    
    methods
        
        %% Constructor
        function obj = sourceorama(srcs,sequenceTimeLength,hdf5file,tel,dataSetTimeLength)
            obj.srcs = srcs;
            obj.sequenceTimeLength = sequenceTimeLength;
            obj.hdf5file = hdf5file;
            obj.tel = tel;
            obj.log = logBook.checkIn(obj);
            if nargin<5
                dataSetTimeLength = 1;
            end
            obj.dataSetTimeLength = dataSetTimeLength;
        end

        %% Destructor
        function delete(obj)
            checkOut(obj.log,obj)
        end
        
        %% Get/Set dataSetTimeLength
        function set.dataSetTimeLength(obj,val)
            obj.p_dataSetTimeLength = val;
            obj.dataSetMemSize = 8*...
                length(obj.srcs)*...
                obj.tel.resolution^2*...
                obj.p_dataSetTimeLength/obj.tel.samplingTime/2^30;
            add(obj.log,obj,sprintf('%3.1fs of data corresponds to %3.1fGB!',obj.p_dataSetTimeLength,obj.dataSetMemSize))
        end
        function val = get.dataSetTimeLength(obj)
            val = obj.p_dataSetTimeLength;
        end
        
        function record(obj)
            
            m_srcs               = obj.srcs;
            m_tel                = obj.tel;
            m_sequenceTimeLength = obj.sequenceTimeLength;
            m_samplingTime       = m_tel.samplingTime;
            m_dataSetTimeLength  = obj.p_dataSetTimeLength;
            m_hdf5file           = obj.hdf5file;
            
            add(obj.log,obj,'Start recording!')

            nDataSet = ceil(m_sequenceTimeLength/m_dataSetTimeLength);
            nSrcs    = length(m_srcs);
            nPx      = m_tel.resolution;
            data     = ...
                zeros(nPx,nPx,nSrcs,ceil(m_dataSetTimeLength/m_samplingTime));
            chronometer = 0;
            m_srcs = m_srcs.*m_tel;
            
            for kDataSet = 1:nDataSet
                
                dataSetTimeLimit = min(kDataSet*m_dataSetTimeLength,m_sequenceTimeLength);
                count = 0;
                
                while chronometer<=dataSetTimeLimit
                    chronometer = chronometer + m_samplingTime;
                    count = count + 1;
                    +m_tel; %#ok<VUNUS>
                    +m_srcs; %#ok<VUNUS>
                    data(:,:,:,count) = cat(3,m_srcs.phase);
                end  
                
                add(obj.log,obj,sprintf('Save batch%d to %s!',kDataSet,m_hdf5file))
                dataSetName = sprintf('/batch%d',kDataSet);
                h5create(m_hdf5file',dataSetName,[nPx,nPx,nSrcs,count]);
                h5write(m_hdf5file',dataSetName,data(:,:,:,1:count));
                
            end
            
            h5writeatt(m_hdf5file,'/','creationDate',datestr(now));
            h5writeatt(m_hdf5file,'/','sequenceTimeLength',m_sequenceTimeLength);
            h5writeatt(m_hdf5file,'/','samplingTime',m_samplingTime);
            
            add(obj.log,obj,'Stop recording!')
        end
        
    end
    
end