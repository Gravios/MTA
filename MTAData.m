classdef MTAData < hgsetget

    properties 
        filename
        path        % path to file containing object
        data        % data
        type        % TimeSeries,TimePeriods,TimePoints,SpatialBins
        ext
        sampleRate  % Sampling rate of data        
    end

    methods
        function Data = MTAData(path,filename,data,sampleRate,type,ext)
            Data.filename = filename;
            Data.path = path;
            Data.data = data;
            Data.type = type;
            Data.ext = ext;
            Data.sampleRate = sampleRate;
        end
        function out = save(Data,varargin)
            out = false;
            data = Data.data;
            if isempty(varargin),
                data = Data.data;
                save(Data.path,'data','-v7.3');
                out = true;
            else
                data = Data.data;
                save(varargin{1},'data','-v7.3');
                out = true;
            end
        end
        function Data = load(Data,varargin)
            if isempty(varargin),
                load(Data.fpath)
                Data.data = data;
            else
                Data.filename = varargin{1};
                load(Data.fpath)
                Data.data = data;
            end
        end
        
        function fpath = fpath(Data)
            fpath = fullfile(Data.path,Data.filename);
        end
        
        function Data = updatePath(Data,path)
            Data.path = path;
        end
        function Data = updateFilename(Data,path)
            Data.path = path;
        end
        
        function sdim = size(Data,varargin)
            if ~isempty(varargin),
                sdim = size(Data.data,cell2mat(varargin));
            else
                sdim = size(Data.data);
            end
        end

        function Data = clear(Data)
            Data.data = [];
        end
        
        function out = isempty(Data)
            out = isempty(Data.data);
        end
        
        function Data = subsref(Data,S)
            ni = numel(S);
            if strcmp(S(1).type,'()')&&ni==1,
                if numel(S.subs)==1,
                    Data = Data.data(S.subs);
                elseif size(S.subs{1},1)==1&&size(S.subs{1},2)>2||...
                       size(S.subs{1},2)==1&&size(S.subs{1},1)>2||...
                       strcmp(S.subs{1},':'),
                    Data = builtin('subsref',Data.data,S);
                else                    
                    Data = SelectPeriods(builtin('subsref',Data.data,S),S.subs{1},'c');
                end
                return
            end
            for n = 1:ni,
                if isa(Data,'MTAData'),
                    dbreak = false;
                    if ismethod(Data,S(n).subs),dbreak = true;end
                    Data = builtin('subsref',Data,S(n:end));
                    if dbreak,break,end
                else
                    Data = subsref(Data,S(n));
                end
            end
        end

    end

    methods (Abstract)        
        Data = create(Data,varargin);      
        Data = resample(data,newSampleRate);
    end

end