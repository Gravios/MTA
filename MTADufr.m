       
classdef MTADufr < MTAData
    properties 
        model = [];
    end
    properties(Transient=true)
        data        % data
    end
    methods
        function Data = MTADufr(varargin)
            [path,filename,data,sampleRate,type,ext] = ...
                DefaultArgs(varargin,{[],[],[],[],'TimeSeries','ufr'});
            if ~isempty(filename),
                if ~strcmp(filename(end-3:end),'.mat'),
                    filename = [filename '.' ext '.mat'];
                end
            end
            Data = Data@MTAData(path,filename,data,sampleRate,type,ext);
        end
        
        function Data = create(Data,Session,varargin)
        %create(Data,Session,DataObj,state,units,twin,overwrite)
        %Calculate the instantaneous firing rate of individual units
        %
            [DataObj,state,units,twin,overwrite] = DefaultArgs(varargin,{Session.lfp,[],[],0.05,0});
            
            if isa(DataObj,'MTAApfs'),
                %nv units
                Data.sampleRate = Session.xyz.sampleRate;
                ind = repmat({zeros(Session.xyz.size(1),1)},1,Session.xyz.size(3));
                for n = 1:numel(DataObj.adata.bins),
                    [~,ind{n}] = NearestNeighbour(DataObj.adata.bins{n}',Session.xyz(:,Session.trackingMarker,n)');
                end
                Data.data = zeros(Session.xyz.size(1),numel(units));
                c = 1;
                for u = units,
                    rateMap = DataObj.data.rateMap(:,ismember(DataObj.data.clu,u),1);
                    Data.data(:,c) = rateMap(sub2ind(DataObj.adata.binSizes',ind{:}));
                    c = c+1;
                end
            
            else
                Data.sampleRate = DataObj.sampleRate;
                syncPeriods = Session.sync.periods(DataObj.sampleRate);
                %if ~exist(Data.fpath,'file')||overwrite,
                spk = Session.spk.copy();
                spk.create(Session,DataObj.sampleRate,state,units);
                Data.data = zeros(diff(syncPeriods([1,end]))+1,numel(units));
                if isempty(units), units = 1:spk.map(end,1);end
                swin = round(twin*DataObj.sampleRate);
                gwin = gausswin(swin)/sum(gausswin(swin));
                for unit = units(:)'
                    Data.data(:,unit==units) = conv(accumarray(spk.res(spk.clu==unit),1,[diff(syncPeriods([1,end]))+1,1]),gwin,'same')/twin;
                end
            end
        end
        
        function Data = load(Data)
        end
        
        function Data = filter(Data)
        end
        function Data = resample(Data,DataObj)
            if DataObj.isempty, DataObj.load; dlen = DataObj.size(1); end
            uind = round(linspace(round(Data.sampleRate/DataObj.sampleRate),Data.size(1),DataObj.size(1)));
            Data.data = Data.data(uind,:);
            Data.sampleRate = DataObj.sampleRate;            
        end
        function Data = embed(Data)
        end

    end
    
end
