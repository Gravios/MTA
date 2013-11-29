classdef MTADepoch < MTAData
    properties
        data 
        label
        key
        model = [];
    end
    
    methods
        
        function Data = MTADepoch(varargin)
            [path,filename,data,sampleRate,label,key,type,ext] = ...
                DefaultArgs(varargin,{[],[],[],[],[],[],'TimePeriods','sst'});
            if ~isempty(filename),
                if ~strcmp(filename(end-3:end),'.mat'),
                    filename = [filename '.' ext '.' label '.' key '.mat'];
                end
            end
            Data = Data@MTAData(path,filename,data,sampleRate,type,ext);
            Data.label = label;
            Data.key = key;
        end
        
        function Data = load(Data,varargin)
            [sync] = DefaultArgs(varargin,{[]});
            if ~isempty(Data.filename)
                load(Data.fpath);
            else
                files = dir(Data.path);
                re = '\.sst\.';
                stsFileList = {files(~cellfun(@isempty,regexp({files.name},re))).name};                
                if ~isempty(Data.key)
                     re = ['\.' Data.key '\.'];
                elseif ~isempty(Data.label)
                    re = ['\.' Data.label '\.'];                   
                else
                    error('STS file does not exist or key/label is missing')
                end
                Data.filename = stsFileList{~cellfun(@isempty,regexp(stsFileList,re))};
                load(Data.fpath)
            end
            if ~isempty(sync),
                sync.resync(Data);
            end
        end
        
        function out = save(Data,varargin)
            [overwrite] = DefaultArgs(varargin,{0});
            out = false;
            
            if ~exist(Data.fpath,'file'),
                save( Data.fpath,'Data','-v7.3');
                out = true;
            elseif exist(Data.fpath,'file')&&overwrite,
                warning(['Overwriting: ' Data.fpath]);
                out = true;
                save( Data.fpath,'Data','-v7.3');
            else
                warning(['File exists: ' Data.fpath, ' - flag the overwrite option  to save']);
            end
            
        end

        function Data = create(Data,varargin)
        end

        function Data = resample(Data,newSampleRate)
            % Needs some more corrections for resampling
            Data.data = round(Data.data/Data.sampleRate*newSampleRate);
            Data.data(Data.data==0)=1;
            Data.sampleRate = newSampleRate;
        end
        
%         function Data = subsref(Data,S)
%             ni = numel(S);
%             if strcmp(S(1).type,'()')&&ni==1,
%                 if numel(S.subs)==0,
%                     Data = Data.data;
%                 elseif numel(S.subs)==1,
%                     if  strcmp(S.subs{1},':')
%                         Data = Data.data(:);
%                     else
%                         Data = Data.data(S.subs);
%                     end
%                 else
%                     Data = Data.data(S.subs{1},S.subs{2});
%                 end
%                 return
%             end
%             for n = 1:ni,
%                 if isa(Data,'MTAData'),
%                     dreturn = false;
%                     if ismethod(Data,S(n).subs),dreturn = true;end
%                     Data = builtin('subsref',Data,S(n:end));
%                     if dreturn,return,end
%                 else
%                     Data = subsref(Data,S(n));                    
%                 end
%             end
%         end

    end
    
    methods (Static)
        function Data = intersect(DataCell)            
            samplingRates = cellfun(@getfield,DataCell,repmat({'sampleRate'},1,numel(DataCell)));
            msr = max(samplingRates);
            if numel(unique(samplingRates))~=1,
                rsi = find(samplingRates~=msr);
                for i = 1:numel(rsi)
                    DataCell{i}.resample(msr);
                end
            end
            newLabel = ['i_' DataCell{1}.key];
            newKey = num2str(randi([0,9],1));
            newData = DataCell{1}.data;
            DataCell(1) = [];
            while ~isempty(DataCell),
                newLabel = [newLabel DataCell{1}.key];
                newData = IntersectRanges(newData,DataCell{1}.data);
                DataCell(1) = [];
            end            
            
            Data = MTADepoch([],[],newData,msr,newLabel,newKey);
        end
        function join(DataCell)
        end
    end
end