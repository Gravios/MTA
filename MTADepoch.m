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
            
            if ~isempty(filename),
                if ~strcmp(filename(end-3:end),'.mat'),
                    filename = [filename '.' ext '.' label '.' key '.mat'];
                end
            end
            Data = Data@MTAData(path,filename,data,sampleRate, ...
                                syncPeriods,syncOrigin,type,ext,name,label,key);
           if isa(path,'MTADepoch'),
                Data = path;
                return
            end

        end
        
        function Data = load(Data,varargin)
        % function Data = load(Data,varargin)
        % 
        % Loads an MTADepoch object from the location in its path property with the
        % filename found in its filename property
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
        % function out = save(Data,varargin)
        % 
        % saves an MTADepoch object under the location in its path property under the
        % filename found in its filename property
        % 
        % Output:
        %     out: boolean - true if save was successful
        %
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
        % function Data = create(Data,Session,funcHandle,varargin)            
        % 
        % Calls one of the heuristic fuctions for state segmentation
        % 
        % Note: This fuction is not ready for general use
        %
            [data,sampleRate,l,k] = feval(funcHandle,Session);
            [label,key] = DefaultArgs(varargin,{l,k});
            Data.path = Session.spath;
            Data.filename = Session.filebase;
            Data.data = data;
            Data.sampleRate = sampleRate;
            Data.label = label;
            Data.key = key;
        end

        function Data = plus(a,b)
        %function Data = plus(a,b)
        %
        % Inuput takes either a two element vector which indicates a time
        %  shift to the start and stop point of each period within the
        %  MTADepoch object.
        %
        % Note: if addition or subtraction results in a period where the
        %       start precedes the stop then that period will be removed
        %
        % Note: if addition or subtraction results two periods overlaping
        %       they will be merged into a single period
        %
        % Future: Addition of two MTADepochs will return their union
        %       
        % Future: suport TimeSeries versions
        
        if isa(a,'MTADepoch')&&~isa(b,'MTADepoch')
            if strcmp(a.type,'TimePeriods'),
                Data = a.copy;
                b = b*a.sampleRate;
                if prod(size(b) == Data.size),
                    Data.data = Data.data+b;
                else
                    Data.data = bsxfun(@plus,Data.data,b);
                end
                Data.clean();
            elseif strcmp(a.type,'TimeSeries'),
            end
            
        elseif isa(b,'MTADepoch')&&~isa(a,'MTADepoch')
            if strcmp(a.type,'TimePeriods'),
                Data = b.copy;
                a = a*b.sampleRate;
                if prod(size(a) == Data.size),
                    Data.data = Data.data+a;
                else
                    Data.data = bsxfun(@plus,Data.data,a);
                end
                
                Data.clean();
            elseif strcmp(a.type,'TimeSeries'),
            end
            
        elseif isa(a,'MTADepoch')&&isa(b,'MTADepoch')
            if strcmp(a.type,'TimePeriods')&&strcmp(b.type,'TimePeriods'),
                Data = a.copy;
                Data.path = [];
                Data.filename = [];
                Data.label = [a.label '+' b.label];
                Data.key = [];
                if b.sampleRate ~= a.sampleRate, b.resample(a); end
                Data.data = JoinRanges(a.data,b.data);
                
                Data.clean();
            elseif strcmp(a.type,'TimeSeries'),
            end
            
        end
        
        % Join the new periods if they overlap
        if sum((Data.data(2:end,1)-Data.data(1:end-1,2))<=0)>0;
            ndata = Data.data(1,:);
            for i = 2:Data.size(1),
                ndata = JoinRanges(ndata,Data.data(i,:));
            end
            Data.data = ndata;
        end
        if Data.sampleRate~=1,
            Data.data = round(Data.data);
        end
        
        
        end

        function Data = minus(a,b)
        % function Data = minus(a,b)
        % 
        % Add assertion that when two MTADepochs are inputed that
        % they have the same origin
              
            if isa(a,'MTADepoch')&&ismatrix(b)&&~isa(b,'MTADepoch'),
                Data = a.copy;
                if all(size(b)==1),
                    Data.data = bsxfun(@minus,Data.data,b);
                else
                    Data.data = SubstractRanges(Data.data,b);
                end
            elseif isa(b,'MTADepoch')&&ismatrix(a)&&~isa(a,'MTADepoch'),
                Data = b.copy;
                if all(size(a)==1),
                    Data.data = bsxfun(@minus,Data.data,a);
                else
                    Data.data = SubstractRanges(a,Data.data);
                end
            elseif isa(b,'MTADepoch')&&isa(a,'MTADepoch')
                b.resample(a.sampleRate);
                Data = a.copy;
                Data.data = SubstractRanges(Data.data,b.data);
                Data.label = [Data.label '-' b.label];
                Data.key = '';
            end
             
        end

         

        function Data = and(a,b)
        %function Data = and(a,b)
        %
        %If one of the imputs is a normal array then
        %the array elements are assumed to be in the 
        %sampling rate of the MTADepoch object.
        %
        % TODO - add MTADepoch/MTADepoch comparisions
            Data = [];
            if isa(a,'MTADepoch')&&ismatrix(b)
                Data = a.copy;
                switch a.type
                  case 'TimePeriods'
                    Data.data = IntersectRanges(Data.data,b);
                    perDur = diff(Data.data,1,2);
                    Data.data(perDur<=0,:) = [];
                  case 'TimeSeries'
                    Data.data = Data.data&b;
                end
                
            elseif isa(b,'MTADepoch')&&ismatrix(a)
                Data = b.copy;
                switch b.type
                  case 'TimePeriods'
                    Data.data = IntersectRanges(Data.data,a);
                    perDur = diff(Data.data,1,2);
                    Data.data(perDur<=0,:) = [];
                  case 'TimeSeries'
                    Data.data = Data.data&a;
                end

            elseif isa(b,'MTADepoch')&&isa(a,'MTADepoch')
                %Data = 
            end
        end
        
        function Data = clean(Data)
        % function Data = clean(Data)
        % Internal utility to clean up periods after set operations
        %
        % Drop periods with zero or negative duration. eg([5,3] or [[10,10])
        % Drop periods which exist before or conains the origin
        % Drop periods which exceed the end of the sync
        % Truncate periods which terminate after end of the sync
        %
            % Drop periods with zero or negative duration. eg([5,3] or [[10,10])
            perDur = diff(Data.data,1,2);
            Data.data(perDur<=0) = [];
            % Drop periods which exist before or conains the origin
            Data.data(Data.data(:,1)<=0|Data.data(:,2)<=0,:) = [];
            % Drop periods which exceed the end of the sync
            Data.data(Data.data(:,1)>round(Data.sync.data(end)*Data.sampleRate)) = [];
            % Truncate periods which terminate after end of the sync
            Data.data( Data.data(:,1)<round(Data.sync.data(end)*Data.sampleRate)...
                &Data.data(:,2)>round(Data.sync.data(end)*Data.sampleRate),2)...
                = round(Data.sync.data(end)*Data.sampleRate);
        end
        
        function Data = fillgaps(Data,varargin)
        % function Data = fillgaps(Data,gap_size)
        % fills gaps between epochs which are smaller than the specified
        % size
        %
        % varargin:
        %   gap_size: numeric, gap size in seconds
        %
        % NOTE: Only functions for TimePeriods
        %
            if ~isempty(varargin),
                gap_size = varargin{1};
            else
                gap_size = round(.1*Data.sampleRate);
            end
        
        
            switch Data.type
                case 'TimePeriods'
                    if size(Data.data,1)>1,
                        interPerDur = Data.data(2:end,1)-Data.data(1:end-1,2);
                        c = 1;
                        while ~isempty(interPerDur)
                            if interPerDur(1)<gap_size,
                                Data.data(c,:)   = [Data.data(c,1),Data.data(c+1,2)];
                                Data.data(c+1,:) = [];
                            else
                                c = c+1;
                            end
                            interPerDur(1) = [];
                        end
                    end
                case 'TimeSeries'
                    perStart = find(diff(double(Data.data))== 1)+1;
                    perStop =  find(diff(double(Data.data))== -1);
                    
                    if ~isempty(perStart)&&~isempty(perStop),
                        if perStop(1)<perStart(1),perStart = [1;perStart]; end
                        if perStop(end)<perStart(end),perStop = [perStop;size(Data,1)]; end
                        
                        interPerDur = perStart(2:end)-perStop(1:end-1);
                        c = 1;
                        while ~isempty(interPerDur)
                            if interPerDur(1)<gap_size,
                                Data.data(perStop(c):perStart(c+1)) = 1;
                                Data.data(c+1,:) = [];
                            else
                                c = c+1;
                            end
                            interPerDur(1) = [];
                        end
                    end
                otherwise
                    error('MTA:MTADepoch:WhatDidYouDO!!!')

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
            newKey = num2str(randi([0,9],1));
            newData = DataCell{1}.data;
            DataCell(1) = [];
            while ~isempty(DataCell),
                newLabel = [newLabel DataCell{1}.key];
                newData = IntersectRanges(newData,DataCell{1}.data);
                DataCell(1) = [];
            end            
            
            Data = MTADepoch([],[],newData,msr,sync,origin,[],[],[],newLabel,newKey);
        end

        function join(DataCell)
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
            newKey = num2str(randi([0,9],1));
            newData = DataCell{1}.data;
            DataCell(1) = [];

            Data = MTADepoch([],[],newData,msr,sync,origin,[],[],[],newLabel,newKey);
        end

% $$$         function Data = plus(a,b)
% $$$             if isa(a,'MTADepoch')&&isvector(b)
% $$$                 Data = a.copy;
% $$$                 b = b*a.sampleRate;
% $$$                 Data.data = Data.data+repmat(b,[Data.size(1),1]);
% $$$                 perDur = dif(Data.data,1,2);
% $$$                 Data.data(perDur<=0) = [];
% $$$             elseif isa(b,'MTADepoch')&&isvector(a)
% $$$                 Data = b.copy;
% $$$                 a = a*b.sampleRate;
% $$$                 Data.data = Data.data+repmat(a,[Data.size(1),1]);
% $$$                 perDur = dif(Data.data,1,2);
% $$$                 Data.data(perDur<=0) = [];
% $$$             elseif isa(b,'MTADepoch')&&isa(a,'MTADepoch')
% $$$                 %Data = 
% $$$             end
% $$$         end


    end
end
