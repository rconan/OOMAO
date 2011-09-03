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
        % data buffer
        buffer
        opdCount;
        bufferCount;
        log;
    end
    
    methods
        
        %% Constructor
        function obj = sourceorama(hdf5file,srcs,sequenceTimeLength,tel,dataSetTimeLength)
            
            obj.srcs = srcs(:)';
            obj.hdf5file = hdf5file;
            
            if nargin>2
                
                obj.sequenceTimeLength = sequenceTimeLength;
                obj.tel = tel;
                obj.log = logBook.checkIn(obj);
                if nargin<5
                    dataSetTimeLength = 1;
                end
                obj.dataSetTimeLength = dataSetTimeLength;
            else
                
                if exist(obj.hdf5file,'file')>0
                    obj.sequenceTimeLength = h5readatt(obj.hdf5file,'/','sequenceTimeLength');
                    obj.buffer      = h5read(obj.hdf5file,'/opdSet1');
                    obj.opdCount    = 1;
                    obj.bufferCount = 0;
                    obj.log = logBook.checkIn(obj);
              else
                    error('oomao:sourceorama','File %s not found!',bj.hdf5file)
                end
                
            end
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
            %% RECORD Record a source time sequence
            %
            % record(obj)
            
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
                
                tic
                while chronometer<=dataSetTimeLimit
                    chronometer = chronometer + m_samplingTime;
                    count = count + 1;
                    +m_tel; %#ok<VUNUS>
                    +m_srcs; %#ok<VUNUS>
                    data(:,:,:,count) = cat(3,m_srcs.opd);
                end  
                elapsedTime = toc;
                
                add(obj.log,obj,sprintf('Elapsed time: %.2fs: save %d runs to opd%d in %s!',...
                    elapsedTime,count,kDataSet,m_hdf5file))
                dataSetName = sprintf('/opdSet%d',kDataSet);
                h5create(m_hdf5file',dataSetName,[nPx,nPx,nSrcs,count]);
                h5write(m_hdf5file',dataSetName,data(:,:,:,1:count));
                
            end
            
            h5writeatt(m_hdf5file,'/','creationDate',datestr(now));
            h5writeatt(m_hdf5file,'/','sequenceTimeLength',m_sequenceTimeLength);
            h5writeatt(m_hdf5file,'/','samplingTime',m_samplingTime);
            h5writeatt(m_hdf5file,'/','sourceZenith',[m_srcs.zenith]);
            h5writeatt(m_hdf5file,'/','sourceAzimuth',[m_srcs.azimuth]);
            h5writeatt(m_hdf5file,'/','sourceHeight',[m_srcs.height]);
            
            add(obj.log,obj,'Stop recording!')
        end

        function uplus(obj)
            %% UPLUS Play a source time sequence
            %
            % +obj
            
            obj.bufferCount = obj.bufferCount + 1;
            for kSrcs = 1:length(obj.srcs)
                obj.srcs(kSrcs).resetPhase = obj.buffer(:,:,kSrcs,obj.bufferCount);
            end
            if obj.bufferCount==size(obj.buffer,4)
                obj.opdCount    = obj.opdCount + 1;
                add(obj.log,obj,sprintf('Loading /opdSet%d',obj.opdCount));
                obj.buffer      = h5read(obj.hdf5file,sprintf('/opdSet%d',obj.opdCount));
                obj.bufferCount = 0;                
            end
        end
        
    end
    
end