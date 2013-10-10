classdef MTASpk < hgsetget
    properties
        
        data
        
        map
        
        sampleRate
        
        type
        
    end
    methods
        
        %s.spk(state,clu,sampleRate,squeeze)
        %s.spk(s.stc{'rear'},
        
        function Spk = MTASpk(varargin)
            [Spk.data,Spk.map,Spk.sampleRate,Spk.type] = ...
                DefaultArgs(varargin,{[],[],[],'TimePoints'});            
        end
        function Spk = create(Spk,Session,varargin)
            %% Load and resample Res
            [Spk.sampleRate,units,states] = DefaultArgs(varargin,{1,[],{}});
            
            [Res Clu Map] = LoadCluRes(fullfile(Session.spath, Session.name));
            Spk.map = Map;
            Res = Res/Session.sampleRate*Spk.sampleRate;
            if spk.sampleRate~=1
                Res = ceil(Res);
            end
            [Res ind] = SelectPeriods(Res,Session.sync([1,end]),'d',1,1);
            Clu = Clu(ind);
            
            Spk.data = cat(2,Res,Clu);
            
            %% Select Specific units
            if isempty(units),
                cind = true(numel(Res),1);
            else
                cind = find(ismember(Clu,units));
            end            
            Spk.data = Spk.data(cind,:);
            
            if ~isempty(states);
                [~,sind] = SelectPeriods(Spk.data(:,1),Session.stc{states,Spk.sampleRate},'d',1,0);
                Spk.data = Spk.data(sind,:);
            end
        end
        
        function resample(Spk)
        end
        
        
    end
end