classdef MTAData < hgsetget
%MTAData(varargin)
%MTAData(path,filename,data,sampleRate,type,ext)
%
%  MTAData is a superclass of most MTA data types. It is a container for 
%  general data types, which alows dynamic refrencing and the use of a set
%  of common functions.
%
%  Current Data Types: 
%    TimeSeries, TimePeriods, TimePoints, SpatialMap
%
%  Current Subclasses:
%    MTADang, MTADepoch, MTADlfp, MTADufr, MTADxyz
%
%  varargin:
%    
%    path:       string, the directory where the object's data is stored
%
%    filename:   string, the file name of the .mat file which contains the 
%                        objects data
%
%    type:       string, a short string which denotes the type of data held
%                        by the object
%
%    ext:        string, a short, unique string which will be the primary
%                        file identifier
%
%    sampleRate: double, the sampling rate of the associated data
%
%


    properties 

        %path - string: the directory where the object's data is stored
        path        

        %filename - string: the file name of the .mat file which contains the objects data
        filename

        %type - string: a short string which denotes the type of data held by the object
        type        
        
        %ext - string: a short, unique string which will be the primary file identifier
        ext

        %sampleRate - double: Sampling rate of the associated data
        sampleRate  
        
        %syncPeriods - numericArray(period,Start/Stop): Time in seconds indicating where the data fits in the Session
        sync
        
        %syncOrigin - double: time of data origin in seconds with respect to the syncPeriods
        origin = 0;
                
    end
    
    properties( Transient=true )
        treatmentRecord = {};
    end
    
    properties (Abstract)
        %data - numericArray: (TimeIndex,Feature,Dimension)
        data        
        
        %model - MTAModel: holds the names of the features for indexing purposes
        model
    end
    
    methods
        function Data = MTAData(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext)
            Data.filename = filename;
            Data.path = path;
            Data.data = data;
            Data.type = type;
            Data.ext = ext;
            Data.sampleRate = sampleRate;
            Data.sync = syncPeriods;
            Data.origin = syncOrigin;
        end
        function out = save(Data,varargin)        
            out = false;
            data = Data.data;
            if isempty(varargin),
                data = Data.data;
                save(Data.fpath,'data','-v7.3');
                out = true;
            else
                data = Data.data;
                save(varargin{1},'data','-v7.3');
                out = true;
            end
        end
        function Data = load(Data,varargin)
            [Session,syncshift] = DefaultArgs(varargin,{[],0});
            ds = load(Data.fpath);
            switch class(Session)
                case 'MTASession'
                    Data.data = ds.data;
                    Data.sync.sync = Session.sync.copy;
                    Session.resync(Data);                    
                case 'MTATrial'
                    Data.data = ds.data;
                    Data.sync.sync = Session.sync.copy;
                    Session.resync(Data);                    
                case 'double'
                    if ~isempty(Session),
                        mf = matfile(Data.fpath);
                        d = ones(1,5);
                        if Data.isempty,
                            dsize(1:numel(size(mf,'data'))) = size(mf,'data');
                        else
                            dsize(1:numel(Data.size)) = Data.size;
                        end
                        d(1:numel(dsize)) = dsize;
                        for i = 1:size(Session,1),
                            Data.data(Session(i,1):Session(i:2),:,:,:,:) = ds.data((Session(i,1):Session(i:2))-syncshift,1:d(2),1:d(3),1:d(4),1:d(5));
                        end
                    else
                        Data.data = ds.data;
                    end
            end
        end
        
        function fpath = fpath(Data)
        %fpath = fpath(Data)
        % Concatenate the path and filename fields to create the full path
        % to the file containing the object's data
        %
            fpath = fullfile(Data.path,Data.filename);
        end
        function Data = updatePath(Data,path)
        %Data = updatePath(Data,path)        
        %
        % Inputs:
        %   path - string: new path to the directory where the Objects data
        %                  should be stored
        %
        % Outputs:
        %   Data - MTAData: Original object passed to this function with an
        %                   updated path
        %
            Data.path = path;
        end
        function Data = updateFilename(Data,filename)
        %Data = updateFilename(Data,filename)        
        %
        % Inputs:
        %   filename - string: new filename for the passed object
        %
        % Outputs:
        %   Data - MTAData: Original object passed to this function with an
        %                   updated filename
        %
            Data.filename = filename;
        end
        function sdim = size(Data,varargin)
        %sdim = size(Data,dimInd)
        % Returns the size of the MTAData object's "data" field
        %
        % Inputs:
        %   dimInd - numericArray: indicies of the dimensions. If v 
        %                            empty all dimesion sizes are
        %                            returned.
        %
        % Outputs:
        %   sdim - numericArray: contains the size of each dimension
        %
            if ~isempty(cell2mat(varargin)),
                dim = cell2mat(varargin);
                ndim = numel(dim);
                if ndim>1,
                    sdim = zeros(1,ndim);
                    for i = 1:ndim,
                        sdim(i) = size(Data.data,dim(i));
                    end
                else
                    sdim = size(Data.data,dim);
                end
            else
                sdim = size(Data.data);
            end
        end
        function Data = clear(Data)
        %Data = clear(Data)
        % Clear an MTAData object's data field.
            Data.data = [];
            Data.origin = [];
        end
        function out = isempty(Data)
            out = isempty(Data.data);
        end
        function Data = subsref(Data,S)
            ni = numel(S);
            switch Data.type
                case 'TimeSeries',
                    if strcmp(S(1).type,'()')&&ni==1,
                        modelRef = (cellfun(@numel,S(1).subs(2:end))>1).*cellfun(@ischar,S(1).subs(2:end))|cellfun(@iscell,S(1).subs(2:end));
                        if sum(modelRef)>0&&~isempty(Data.model),
                            mri = find(modelRef)+1;
                            for  i = mri,
                                S(1).subs{i} = Data.model.gmi(S(1).subs{i});
                            end
                        end
                        
                        % Convert MTAepoch to nx2 matrix and resample if necessary
                        epicDataInd = find(cellfun(@isa,S(1).subs,repmat({'MTADepoch'},1,numel(S(1).subs))));
                        for i = epicDataInd,
                            epoch = S(1).subs{i}.copy;
                            if Data.sampleRate~=epoch.sampleRate,
                                epoch.resample(Data.sampleRate);
                            end
                            S(1).subs{i} = epoch.data;
                        end
                        
                        if numel(S.subs)==0,
                            Data = Data.data;
                        elseif numel(S.subs)==1,
                            if (size(S.subs{1},1)==1&&size(S.subs{1},2)>2)||...
                                    (size(S.subs{1},2)==1&&size(S.subs{1},1)>2)||...
                                    strcmp(S.subs{1},':'),
                                
                                Data = builtin('subsref',Data.data,S);
                                return
                            else
                                S.subs = S.subs(~cellfun(@isempty,S.subs));
                                Sa = S;
                                Sa.subs{1} = ':';
                                Data = SelectPeriods(builtin('subsref',Data.data,Sa),S.subs{1},'c');
                            end
                        elseif (size(S.subs{1},1)==1&&size(S.subs{1},2)>2)||...
                                (size(S.subs{1},2)==1&&size(S.subs{1},1)>=2)||...
                                strcmp(S.subs{1},':'),
                            Data = builtin('subsref',Data.data,S);
                        else
                            S.subs = S.subs(~cellfun(@isempty,S.subs));
                            Sa = S;
                            Sa.subs{1} = ':';
                            Data = SelectPeriods(builtin('subsref',Data.data,Sa),S.subs{1},'c');
                        end
                        return
                    end
                case 'TimePeriods'
                    if strcmp(S(1).type,'()')&&ni==1,
                        if numel(S.subs)==0,
                            Data = Data.data;
                        elseif numel(S.subs)==1,
                            if  strcmp(S.subs{1},':')
                                Data = Data.data(:);
                            else
                                Data = Data.data(S.subs{1});
                            end
                        else
                            Data = Data.data(S.subs{1},S.subs{2});
                        end
                        return
                    end
            end
            if strcmp(S(1).type,'.'),
                if isprop(Data,S(1).subs)
                    if numel(S)==1,
                        Data = Data.(S(1).subs);
                    else
                        Data = subsref(Data.(S(1).subs),S(2:end));
                    end
                elseif ismethod(Data,S(1).subs)
                    if numel(S)==1,
                        Data = Data.(S(1).subs);
                    else
                        Data = Data.(S(1).subs)(S(2).subs{:});
                    end
                else
                    Data = builtin('subsref',Data,S);
                end
            end

        end
        function DataCopy = copy(Data)
        % Make a copy of a handle object.
        % Instantiate new object of the same class.
            DataCopy = feval(class(Data),[]);
            % Copy all non-hidden properties.
            p = properties(Data);
            for i = 1:length(p)
                if isa(Data.(p{i}),'MTAData'),
                    DataCopy.(p{i}) = Data.(p{i}).copy;
                else
                    DataCopy.(p{i}) = Data.(p{i});
                end
            end
        end
        
        function Data = cast(Data,type,varargin)
            [sampleRate,syncflag] = DefaultArgs(varargin,{Data.sampleRate,'relative'});
            oldType = Data.type;            
            if ~strcmp(oldType,type)&&~isempty(oldType)
                switch type
                    case 'TimePeriods'
                        Data.data = ThreshCross(Data.data,0.5,1);
                    case 'TimeSeries'
                        tmpdata = round(Data.data./Data.sampleRate.*sampleRate)+Data.origin;
                        tmpdata(tmpdata==0) = 1;
                        if isa(Data.sync,'MTAData'),
                            Data.sync.resample(Data.sampleRate);
                            data = false(round(Data.sync.data(end)./Data.sampleRate.*sampleRate),1);
                        else
                            data = false(round(Data.sync(end)./Data.sampleRate.*sampleRate),1);
                        end                                
                        for j = 1:size(tmpdata,1),
                            data(tmpdata(j,1):tmpdata(j,2)) = true;
                        end         
                        switch syncflag
                            case 'relative'
                                try
                                    data = [data(Data.sync.sync.data(1)+1:(Data.sync.sync.data(end)));0];
                                catch err
                                    data = cat(1,data,zeros(Data.sync.sync.data(end)-Data.sync.data(end),1));
                                    data = data(Data.sync.sync.data(1)+1:Data.sync.sync.data(end));
                                end
                            case 'absolute'

                        end
                        
                        Data.data = data;
                end
            end            
            Data.type = type;
        end
        
    end

    methods (Abstract)        
        Data = create(Data,varargin);      
        Data = resample(Data,newSampleRate);
        %Data = filter(Data,fwindow);
    end

end