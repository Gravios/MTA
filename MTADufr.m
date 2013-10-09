       
classdef MTADang < MTAData

    properties(Transient=true)
        data        % data
    end
    methods
        function Data = MTADang(varargin)
            [path,filename,data,sampleRate,type,ext] = ...
                DefaultArgs(varargin,{[],[],[],[],'TimeSeries','ufr'});
            if ~isempty(filename),
                if ~strcmp(filename(end-3:end),'.mat'),
                    filename = [filename '.' ext '.mat'];
                end
            end
            Data = Data@MTAData(path,filename,data,sampleRate,type,ext);
        end
        
        function Data = create(Data,Session)
            [newSampleRate,flength,twin] = DefaultArgs(varargin,{Session.lpf.sampleRate,Session.xyzPeriods(end)-Session.xyzPeriods(1)+1,0.05});
            if ~exist([Session.spath Session.name '.ufr'],'file'),
                [Res,Clu,Map] = LoadCluRes(fullfile(Session.spath,Session.name));

                gwin = gausswin(round(twin*Session.sampleRate))/sum(gausswin(round(twin*Session.sampleRate)));
                spks = zeros(max(Res),1);
                ufr = [];
                for unit = 1:size(Map,1),
                    uRes = Res(Clu==unit);
                    uClu = Clu(Clu==unit);
                    spks(:)=0;
                    spks(uRes) = 1;
                    ufr(:,unit) =resample(conv(spks,gwin),Session.lpf.sampleRate,Session.sampleRate);
                end
                
                save([Session.spath Session.name '.ufr'],'ufr','-v7.3','-mat');
                ufr = ufr(Session.syncPeriods(1):Session.syncPeriods(end),:);
            else 
                load([Session.spath Session.name '.ufr'],'-mat')
                ufr = ufr(Session.syncPeriods(1):Session.syncPeriods(end),:);
            end
            
            if newSampleRate == Session.lpf.sampleRate,
                Session.ufr = ufr;
            else
                uind = round(linspace(round(Session.lpf.sampleRate/newSampleRate),size(ufr,1),flength));
                Session.ufr = ufr(uind,:);
            end
        end
        
        function Data = filter(Data,win)
        end
        function Data = resample(Data,newSampleRate,varargin)
            [interp_type] = DefaultArgs(varargin,{'linear'});
        end
        function Data = embed(Data,win,overlap)
        end

    end
    
end
