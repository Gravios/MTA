classdef MTADufr < MTAData
    properties 
        model = [];
    end
    properties(Transient=true)
        data        % data
    end
    methods
        function Data = MTADufr(varargin)
            [path,filename,data,sampleRate,syncPeriods,syncOrigin,model,type,ext,name,label,key] = ...
                DefaultArgs(varargin,{[],[],[],[],[],[],[],'TimeSeries','ufr',[],'ufr','u'});
            Data = Data@MTAData(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext,name,label,key);
        end
        
        function Data = create(Data,Session,varargin)
        %create(Data,Session,DataObj,state,units,twin,overwrite)
        %Calculate the instantaneous firing rate of individual units
        %
            [DataObj,state,units,twin,overwrite] = DefaultArgs(varargin,{Session.lfp,[],[],0.05,0});
            
            if isa(DataObj,'MTAApfs'),
                %nv units
                if isempty(units),
                    units = cat(1,DataObj.data.clu);
                end
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
                %if ~exist(Data.fpath,'file')||overwrite,
                spk = Session.spk.copy();
                spk.create(Session,DataObj.sampleRate,state,units);
                if ~DataObj.isempty,
                    dsize = DataObj.size(1);
                else
                    dsize = diff(round(DataObj.sync.sync([1,end])*DataObj.sampleRate+1));
                end
                Data.data = zeros(dsize,numel(units));
                if isempty(units), units = 1:spk.map(end,1);end
                swin = round(twin*DataObj.sampleRate);
                %gwin = gausswin(swin)/sum(gausswin(swin));
                gwin = ones(swin,1);
                for unit = units(:)'
                    res = spk.res(spk.clu==unit);
                    res(res==0) = 1;
                    res(res>dsize) = dsize;
                    Data.data(:,unit==units) = conv(accumarray(res,1,[dsize,1])./twin,gwin,'same');
                end
                Data.origin = DataObj.origin;
                Data.sync = DataObj.sync;
            end
        end
        
        function Data = load(Data)
        end
        
        function Data = embed(Data)
        end

    end
    
end
