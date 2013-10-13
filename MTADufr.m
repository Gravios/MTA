       
classdef MTADufr < MTAData

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
            
            [DataObj,twin] = DefaultArgs(varargin,{Session.lfp,0.05});
           
            lfpSyncPeriods = Session.sync.periods(Session.lfp.sampleRate);
            
            if ~exist(Data.fpath,'file'),
                [Res,Clu,Map] = LoadCluRes(fullfile(Session.spath,Session.name));

                gwin = gausswin(round(twin*Session.sampleRate))/sum(gausswin(round(twin*Session.sampleRate)));
                spks = zeros(max(Res),1);
                ufr = [];
                for unit = 1:size(Map,1),
                    uRes = Res(Clu==unit);
                    uClu = Clu(Clu==unit);
                    spks(:)=0;
                    spks(uRes) = 1;
                    ufr(:,unit) = resample(conv(spks,gwin),Session.lfp.sampleRate,Session.sampleRate);
                end
                
                save(Data.fpath,'ufr','-v7.3','-mat');
            else 
                load([Session.spath Session.name '.ufr'],'-mat')
            end
            
            ufr = ufr(lfpSyncPeriods(1):lfpSyncPeriods(end),:);
            
            if DataObj.sampleRate < Session.lpf.sampleRate,
                Data.resample(DataObj);
            end

        end
        
        function Data = load(Data)
        end
        
        function Data = filter(Data)
        end
        function Data = resample(Data,DataObj)
            if DataObj.isepmty, DataObj.load; dlen = DataObj.size(1); end
            uind = round(linspace(round(Data.sampleRate/DataObj.sampleRate),Data.size(1),DataObj.size(1)));
            Data.data = Data.data(uind,:);
            Data.sampleRate = DataObj.sampleRate;
            DataObj.clear;
        end
        function Data = embed(Data)
        end

    end
    
end
