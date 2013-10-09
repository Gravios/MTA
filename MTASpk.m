classdef MTASpk < hgsetget
    properties
        
        clu
        
        res
        
        map
        
        sampleRate
        
    end
    methods
        
        function Spk = MTASpk(varargin)
            [res,clu,map,sampleRate] = ...
                DefaultArgs(varargin,{[],[],[],[]});
            
        end
        function Spk = create(Spk,Session,varargin)
            %Session = load_CluRes(Session,varargin)
            %load unit timing
            %sampleRate - double: sample rate to which res is to be resampled
            %  Default: 1250 - lfpSampleRate
            %
            %units - int: the unit id's you wish to load.
            %states - string/cellArray(string): specify a state or a
            %         union/intersection of states
            %  note: only use single states for the moment!!
            %  Default: {}
            %
            %sets:
            %    Session.res = Res;
            %    Session.clu = Clu;
            %    Session.map = Map;            
            
            [sampleRate,units,states] = DefaultArgs(varargin,{Session.lpf.sampleRate,[],{}});
            [Res Clu Map] = LoadCluRes(fullfile(Session.spath, Session.name));
            Res = round(Res*sampleRate/Session.sampleRate);
            [Res ind] = SelectPeriods(Res,round([Session.syncPeriods(1),Session.syncPeriods(end)]./Session.lpf.sampleRate.*sampleRate),'d',1,1);
            Clu = Clu(ind);
            if isempty(units),
                cind = true(numel(Res),1);
            else
                cind = find(ismember(Clu,units));
            end
            Spk.res = Res(cind);
            Spk.clu = Clu(cind);
            Spk.map = Map;
            if ~isempty(states);
                [Spk.res sind] = SelectPeriods(Spk.res,Session.statePeriods(states,sampleRate),'d',1,0);
                Session.clu = Session.clu(sind);
            end
        end
        
    end
end