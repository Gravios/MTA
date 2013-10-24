classdef MTASpk < hgsetget
    properties
        
        clu
        
        res
        
        map
        
        sampleRate
        
        type
        
    end
    methods
        
        %s.spk(state,clu,sampleRate,squeeze)
        %s.spk(s.stc{'rear'},
        
        function Spk = MTASpk(varargin)
            [Spk.res,Spk.clu,Spk.map,Spk.sampleRate,Spk.type] = ...
                DefaultArgs(varargin,{[],[],[],[],'TimePoints'});            
        end
        function Spk = create(Spk,Session,varargin)
            %% Load and resample Res             
            [Spk.sampleRate,states,units] = DefaultArgs(varargin,{1,[],{}});
            
            [Res, Clu, Map] = LoadCluRes(fullfile(Session.spath, Session.name));
            Spk.map = Map;
            Res = Res/Session.sampleRate*Spk.sampleRate;
            if Spk.sampleRate~=1
                Res = ceil(Res);
            end
            [Res, ind] = SelectPeriods(Res,ceil(Session.sync([1,end])*Spk.sampleRate),'d',1,1);
            Clu = Clu(ind);
            
            %% Select specific units
            if isempty(units),
                cind = true(numel(Res),1);
            else
                cind = find(ismember(Clu,units));
            end            
            
            Res = Res(cind);
            Clu = Clu(cind);
            
            %% Select specific states
            if ~isempty(states);
                [Res,sind] = SelectPeriods(Res,[Session.stc{states,Spk.sampleRate}],'d',1,0);
                Clu = Clu(sind);
            end
            
            Spk.clu = Clu;
            Spk.res = Res;
        end
        
        function Data = subsref(Data,S)
            ni = numel(S);
            if strcmp(S(1).type,'()')&&ni==1,
                if numel(S.subs)==0,
                    Data = Data.res;
                else
                    Data = Data.res(ismember(Data.clu,S.subs{1}));
                end
                return
            end
            Data = builtin('subsref',Data,S);
        end
        
        function Data = clear(Data)
            Data.res = [];
            Data.clu = [];
            Data.sampleRate = [];
        end
        
        
    end
end