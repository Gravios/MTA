classdef MTADepoch < MTAData
% MTADepoch(path,filename,data,sampleRate,syncOrigin,label,key,datatype,extension)
% MTADepoch is a container for periods of time as defined by start and stop
% instances recorded in a Nx2 array, within the context of a dynamic
% temporal domain (e.g. MTASession or MTATrial).
%
%  Data Types:
%
%    TimePeriods: Nx2 array, The index boundaries of the stored epochs
%    TimeSeries:  Tx1 array, Logical array where, (1 == within epoch)
%
%  Behavior:
%    
%      
%   
%  Variables:
%    
%    path:        string, location where the data should be saved
%    
%    filename:    string, Name of the file
%
%    data:        double, matrix of size Nx2, where N is the number of epochs
%
%    sampleRate:  double, Number of samples per second to which start and stop
%                         values of the data corespond.
%
%    syncPeriods: MTADepoch, time periods in reference to larger timescale
%                            indicating where the data fits in a Session

    properties
        data 
        label
        key
        model = [];
    end
    
    methods
        
        function Data = MTADepoch(varargin)
            [path,filename,data,sampleRate,syncPeriods,syncOrigin,label,key,type,ext] = ...
                DefaultArgs(varargin,{[],[],[],[],[],[],[],[],'TimePeriods','sst'});
            if ~isempty(filename),
                if ~strcmp(filename(end-3:end),'.mat'),
                    filename = [filename '.' ext '.' label '.' key '.mat'];
                end
            end
            Data = Data@MTAData(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext);
            Data.label = label;
            Data.key = key;
        end
        
        function Data = load(Data,varargin)
            [Session,sync] = DefaultArgs(varargin,{[],[]});
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
            if ~isempty(Session),
                Session.resync(Data);
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

        function Data = create(Data,Session,funcHandle,varargin)            
            [data,sampleRate,l,k] = feval(funcHandle,Session);
            [label,key] = DefaultArgs(varargin,{l,k});
            Data.path = Session.spath;
            Data.filename = Session.filebase;
            Data.data = data;
            Data.sampleRate = sampleRate;
            Data.label = label;
            Data.key = key;
        end

        function Data = resample(Data,newSampleRate)
            % Needs some more corrections for resampling
            if newSampleRate == 1
                rf = @(x)x;
            elseif isa(newSampleRate,'MTAData')
                newSampleRate = newSampleRate.sampleRate;
                rf = @round;
            else
                rf = @round;
            end
            
            Data.data = rf(Data.data/Data.sampleRate*newSampleRate);
            while sum(Data.data(:)==0)>1
                Data.data(1,:) = [];
            end
            if isa(Data.sync,'MTAData'),
                Data.sync.resample(newSampleRate);
                Data.origin = rf(Data.origin/Data.sampleRate*newSampleRate); 
            else
                Data.sync = rf(Data.sync/Data.sampleRate*newSampleRate);
                Data.sync(Data.sync==0) = 1;
            end            
            Data.origin(Data.origin==0) = 1;
            Data.data(Data.data==0)=1;
            Data.sampleRate = newSampleRate;
        end
        
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