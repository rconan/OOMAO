classdef telescoporama < handle
    %% SOURCEORAMA Create a sourceorama object
    %
    % srcorma =
    % sourceorama(hdf5file,srcs,sequenceTimeLength,tel,dataSetTimeLength)
    % create a sourceorama object from a hdf5 filename, a sources vector,
    % the time length of the source propagation sequence, a
    % telescope+atmosphere object and the time length of one data set
    %
    % srcorma = sourceorama(hdf5file,srcs)
    % create a sourceorama object from a hdf5 filename and a sources
    % vector. This constructor is used to re-play a recorded source
    % propagation sequence.
   
    properties
        
        % source object array
        srcs;
        % atmosphere turbulence layers
        layers;
        % telescope+atmosphere object
        tel;
        % HDF5 file name
        hdf5file;
        % whole data time length
        sequenceTimeLength
        % tag
        tag = 'TELESCOPORAMA';
        newh5File = true;
        chronometer = 0;
    end
    
    properties (Dependent)
        % time length of each data set [s]
        dataSetTimeLength;
        
    end
    
    properties (SetAccess=private)
        % memory size of each data set [GB]
        dataSetMemSize = 0.1;
    end
    
    properties (Access=private)
        % time length of each data set [s]
        p_dataSetTimeLength;
        p_dataSetLength;
        % data buffer
        buffer
        opdCount = 1;
        bufferCount = 0;;
        dataSetTimeLimit;
        log;
    end
    
    methods
        
        %% Constructor
        function obj = telescoporama(hdf5file,tel,varargin)
            
            p = inputParser;
            p.addRequired('hdf5file'          , @ischar);
            p.addRequired('tel'               , @(x) isa(x,'telescopeAbstract') );
            p.addOptional('sequenceTimeLength', Inf, @isnumeric );
            p.addOptional('dataSetTimeLength' , 1  , @isnumeric );
            p.parse(hdf5file,tel,varargin{:});
            
            obj.hdf5file = p.Results.hdf5file;
            obj.sequenceTimeLength ...
                = p.Results.sequenceTimeLength;
            obj.tel      = p.Results.tel;
            obj.layers   = tel.opticalAberration.layer; 
            % May be needed by does not seem to work!
%             obj.tel      = obj.tel - obj.tel.opticalAberration;
%             obj.srcs     = obj.srcs.*obj.tel;
%             obj.tel      = obj.tel + obj.tel.opticalAberration;
            obj.log      = logBook.checkIn(obj);
            obj.dataSetTimeLength ...
                = p.Results.dataSetTimeLength;

            if exist(obj.hdf5file,'file')>0
                add(obj.log,obj,sprintf('Found HDF5 file %s, data will be extracted from it!',...
                    obj.hdf5file))
                disp(h5info(obj.hdf5file))
                obj.newh5File = false;
                obj.sequenceTimeLength = h5readatt(obj.hdf5file,'/','sequenceTimeLength');
                obj.buffer      = h5read(obj.hdf5file,'/opdSet1');
                obj.log = logBook.checkIn(obj);
            end
            
        end

        %% Destructor
        function delete(obj)
            if obj.newh5File
                
                m_tel                = obj.tel;
                m_hdf5file           = obj.hdf5file;
                
                if obj.chronometer <= obj.sequenceTimeLength
                    nPx      = m_tel.resolution;
                    add(obj.log,obj,sprintf('Save %d runs to opd%d in %s!',...
                        obj.bufferCount,obj.opdCount,m_hdf5file))
                    for k=1:length(obj.buffer)
                        dataSetName = sprintf('/opdSet%d/layer%d',obj.opdCount,k);
                        h5create(m_hdf5file',dataSetName,[nPx,nPx,obj.bufferCount]);
                        h5write(m_hdf5file',dataSetName,obj.buffer{k}(:,:,1:obj.bufferCount));
                    end
                end
                
                h5writeatt(m_hdf5file,'/','creationDate',datestr(now));
                h5writeatt(m_hdf5file,'/','sequenceTimeLength',obj.sequenceTimeLength);
                h5writeatt(m_hdf5file,'/','samplingTime',m_tel.samplingTime);
            end
            checkOut(obj.log,obj)
        end
        
        %% Get/Set dataSetTimeLength
        function set.dataSetTimeLength(obj,val)
            obj.dataSetMemSize = 8*...
                length(obj.srcs)*...
                obj.tel.resolution^2*...
                obj.p_dataSetTimeLength/obj.tel.samplingTime/2^30;
            obj.dataSetMemSize = 0.01;
            obj.p_dataSetLength = max( round( obj.dataSetMemSize*2^30/( 8*sum([obj.layers.nPixel].^2) ) ) , 1);
            obj.p_dataSetTimeLength = obj.tel.samplingTime*obj.p_dataSetLength;
            add(obj.log,obj,sprintf('%.3fs of data corresponds to %.2fGB!',obj.p_dataSetTimeLength,obj.dataSetMemSize))
            
            obj.dataSetTimeLimit = min(obj.opdCount*obj.p_dataSetTimeLength,obj.sequenceTimeLength);
            
            if ~isempty(obj.tel)
                obj.buffer = cellfun( @(x) zeros(x,x,obj.p_dataSetLength) , {obj.layers.nPixel} , 'UniformOutput', false);
            end
            
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
            data     = obj.buffer;
            obj.buffer = [];
            m_chronometer = 0;
            m_srcs = m_srcs.*m_tel;
            
            for kDataSet = 1:nDataSet
                
                m_dataSetTimeLimit = min(kDataSet*m_dataSetTimeLength,m_sequenceTimeLength);
                count = 0;
                
                tic
                while m_chronometer<m_dataSetTimeLimit
                    m_chronometer = m_chronometer + m_samplingTime;
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
            obj.newh5File = false;
            
        end

        function uplus(obj)
            %% UPLUS Play a source time sequence
            %
            % +obj
            
            if obj.newh5File
                
                if obj.chronometer<=obj.sequenceTimeLength
                    
                    m_tel                = obj.tel;
                    
                    obj.chronometer = obj.chronometer + obj.tel.samplingTime;
                    
                    obj.bufferCount  = obj.bufferCount + 1;
                    +m_tel; %#ok<VUNUS>
                    
                    for k=1:length(obj.buffer)
                        obj.buffer{k}(:,:,obj.bufferCount) = obj.layers(k).phase;
                    end
                    
                    if obj.chronometer>=obj.dataSetTimeLimit
                        
                        m_hdf5file           = obj.hdf5file;
                        
                        add(obj.log,obj,sprintf('Save %d runs to opd%d in %s!',...
                            obj.bufferCount,obj.opdCount,m_hdf5file))
                        for k=1:length(obj.buffer)
                            dataSetName = sprintf('/opdSet%d/layer%d',obj.opdCount,k);
                            nPx      = obj.layers(k).nPixel;;
                            h5create(m_hdf5file',dataSetName,[nPx,nPx,obj.bufferCount]);
                            h5write(m_hdf5file',dataSetName,obj.buffer{k}(:,:,1:obj.bufferCount));
                        end
                        
                        obj.opdCount = obj.opdCount + 1;
                        obj.dataSetTimeLimit = min(obj.opdCount*obj.p_dataSetTimeLength,obj.sequenceTimeLength);
                        obj.bufferCount = 0;
                        
                    end
                    
                end
                
            else
                
                if obj.bufferCount==size(obj.buffer{1},3)
                    obj.opdCount    = obj.opdCount + 1;
                    for k=1:length(obj.buffer)
                        add(obj.log,obj,sprintf('Loading /opdSet%d/layer%d',obj.opdCount,k));
                        obj.buffer{k}      = h5read(obj.hdf5file,sprintf('/opdSet%d/layer%d',obj.opdCount,k));
                    end
                    obj.bufferCount = 0;
                end
                obj.bufferCount = obj.bufferCount + 1;
                for kSrcs = 1:length(obj.srcs)
                    for k=1:length(obj.buffer)
                        obj.layers(k).phase = obj.buffer{k}(:,:,obj.bufferCount);
                    end
                end
                
            end
        end
        
        function reset(obj)
            obj.chronometer = 0;
            obj.bufferCount = 0;
            obj.opdCount    = 1;
        end
        
    end
    
end