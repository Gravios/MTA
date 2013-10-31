       
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
            
            [DataObj,units,twin] = DefaultArgs(varargin,{Session.lfp,[],0.05});
            
            if isa(DataObj,'MTAApfs'),
                %nv units
                ind = repmat({zeros(Session.xyz.size(1),1)},1,Session.xyz.size(3));
                for n = 1:numel(DataObj.adata.bins),
                    [~,ind{n}] = NearestNeighbour(DataObj.adata.bins{n}',Session.xyz(:,Session.trackingMarker,n)');
                end
                
                Data.data = zeros(Session.xyz.size(1),numel(units));
                c = 1;
                for u = units,
                    rateMap = DataObj.data.rateMap(:,ismember(DataObj.data.clu,units),1);
                    Data.data(:,c) = rateMap(sub2ind(DataObj.adata.binSizes',ind{:}));
                    c = c+1;
                end
                Data.sampleRate = Session.xyz.sampleRate;
                
            else
                lfpSyncPeriods = Session.sync.periods(Session.lfp.sampleRate);
                if ~exist(Data.fpath,'file'),
                    [Res,Clu,Map] = LoadCluRes(fullfile(Session.spath,Session.name));
                    
                    gwin = gausswin(round(twin*Session.sampleRate))/sum(gausswin(round(twin*Session.sampleRate)));
                    spks = zeros(max(Res),1);
                    data = [];
                    for unit = 1:size(Map,1),
                        uRes = Res(Clu==unit);
                        uClu = Clu(Clu==unit);
                        spks(:)=0;
                        spks(uRes) = 1;
                        data(:,unit) = resample(conv(spks,gwin),Session.lfp.sampleRate,Session.sampleRate);
                    end
                    
                    save(Data.fpath,'data','-v7.3','-mat');
                else
                    load(Data.fpath);
                end
                
                if isempty(units),units = ':';end
                Data.data = data(lfpSyncPeriods(1):lfpSyncPeriods(end),units);
             
                if DataObj.sampleRate < Session.lfp.sampleRate,
                    Data.resample(DataObj);
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
            DataObj.clear;
        end
        function Data = embed(Data)
        end

    end
    
end
