classdef MTADepoch < MTAData
% MTADepoch(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext,name,label,key)
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
%  Operators:
%
%    "+": MTADepoch + matrix (matrix values must be in seconds)
%           (TimePeriods) Input: 1x2 matrix, which contains an amount of
%                           time to add to either end of each start/stop 
%                           within the collection of stored epochs
%           (TimePeriods) Input: Nx2 matrix, the MTADepoch must have N rows
%                           Then MTADepoch and the Array will be added normally
%
%    "+": MTADepoch + MTADepoch (NOT FUNCTIONAL)
%           (TimePeriods) Future, Will use JoinRanges to form the union of
%                           the two epoch collections
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
%
%    syncOrigin:  double, point of origin in the overall session timeline
%                         where the obejec begins (rewrite this)
%
%    label:       string, name associated with the type of epochs
%
%    key:           char, single character used for keyboard shortcuts and indexing
%

    properties
        %data - double: matrix of size Nx2, where N is the number of epochs
        data 
        
        %model - Reserved for future versions
        model = [];
    end
    
    methods
        
        function Data = MTADepoch(varargin)
            [path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext,name,label,key] = ...
                DefaultArgs(varargin,{[],[],[],[],[],[],'TimePeriods','sst','epochs',[],[]});
            
            Data = Data@MTAData(path,filename,data,sampleRate, ...
                                syncPeriods,syncOrigin,type,ext,name,label,key);
            if isa(path,'MTADepoch'),
                Data = path;
                return
            end

        end
     
        
    end        
    
    methods (Static)
        function Data = intersect(DataCell)            
            samplingRates = cellfun(@getfield,DataCell,repmat({'sampleRate'},1,numel(DataCell)));
            msr = max(samplingRates);
            if numel(unique(samplingRates))~=1,
                rsi = find(samplingRates~=msr);
                for i = rsi,
                    DataCell{i}.resample(msr);
                end
            end
            sync = DataCell{1}.sync.copy;
            origin = DataCell{1}.origin;
            newLabel = ['i_' DataCell{1}.key];
            %newKey = num2str(randi([0,9],1));
            newKey = 'c';
            newData = DataCell{1}.data;
            DataCell(1) = [];
            while ~isempty(DataCell),
                newLabel = [newLabel DataCell{1}.key];
                newData = IntersectRanges(newData,DataCell{1}.data);
                DataCell(1) = [];
            end            
            
            Data = MTADepoch([],[],newData,msr,sync,origin,[],[],[],newLabel,newKey);
        end

        function Data = join(DataCell)
            if numel(DataCell)==1||isa(DataCell,'MTADepoch'),
                if iscell(DataCell)
                    Data = DataCell{1};
                else
                    Data = DataCell;
                end
                return;
            end
            samplingRates = cellfun(@getfield,DataCell,repmat({'sampleRate'},1,numel(DataCell)));
            msr = max(samplingRates);
            if numel(unique(samplingRates))~=1,
                rsi = find(samplingRates~=msr);
                for i = rsi,
                    DataCell{i}.resample(msr);
                end
            end
            for i = 1:numel(DataCell)
                DataCell{i}.cast('TimeSeries');
            end

            sync = DataCell{1}.sync.copy;
            origin = DataCell{1}.origin;

            newLabel = ['u_' DataCell{1}.key];
            %newKey = num2str(randi([0,9],1));
            newKey = 'c';
            newData = DataCell{1}.data;
            DataCell(1) = [];

            Data = MTADepoch([],[],newData,msr,sync,origin,[],[],[],newLabel,newKey);
        end



    end
end
