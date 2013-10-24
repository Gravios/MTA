classdef MTAData < hgsetget

    properties 
        filename
        path        % path to file containing object
        type        % TimeSeries,TimePeriods,TimePoints,SpatialBins
        ext
        sampleRate  % Sampling rate of data
    end
    
    properties (Abstract)
        data        % data
        model
    end
    
    methods
        function Data = MTAData(varargin)
            [path,filename,data,sampleRate,type,ext] = ...
                DefaultArgs(varargin,{[],[],[],[],[],[]});
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
            [periods] = DefaultArgs(varargin,{[]});
            load(Data.fpath)
            if ~isempty(periods),
                Data.data = data(periods(1):periods(2),:,:,:,:);
            else
                Data.data = data;
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
            if ~isempty(varargin)&&isempty(cell2mat(varargin)),
                sdim = size(Data.data,cell2mat(varargin));
            else
                sdim = size(Data.data);
            end
        end

        function Data = clear(Data)
        %Data = clear(Data)
        % Clear an MTAData object's data field.
            Data.data = [];
        end
        
        function out = isempty(Data)
            out = isempty(Data.data);
        end
        
        function Data = subsref(Data,S)
            ni = numel(S);
            if strcmp(S(1).type,'()')&&ni==1,
                modelRef = (cellfun(@numel,S(1).subs(2:end))>1).*cellfun(@ischar,S(1).subs(2:end));
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
                   
                        Data = Data.data(S.subs{1});
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
            for n = 1:ni,
                if isa(Data,'MTAData'),
                    Data = builtin('subsref',Data,S(n:end));
                    break
                else
                    Data = subsref(Data,S(n));
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
                DataCopy.(p{i}) = Data.(p{i});
            end
        end
        
    end

    methods (Abstract)        
        Data = create(Data,varargin);      
        Data = resample(data,newSampleRate);
    end

end