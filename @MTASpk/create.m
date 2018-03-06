function Spk = create(Spk,Session,varargin)
% Spk = create(Spk,Session,varargin)
% Function used to load spiking data from {clu,res,fet,spk} file
% and synchronize with the Session
%   
%   varargin:
%     [sampleRate,states,units,mode,loadField] 
%
%     sampleRate: double, the sampleRate in Hertz
%
%     states:     string, state label or expression to select for a
%                         subset of spikes
%
%     units:       array, Cluster identities to select for a subset
%                         of spikes. The default [] will return all
%                         Clusters.
%
%     mode:       string, {'','deburst'}, only keep 1st spk of each
%                         burst.
%
%     loadField: cellarry, NOT IMPLEMENTED
%
% See also MTAStateCollection for more information on selecting
% states
%
%% Load and resample Res             
[sampleRate,states,units,mode,loadField] = DefaultArgs(varargin,{1,[],[],'',{}});
Spk.sampleRate = sampleRate;
[Res, Clu, Map] = LoadCluRes(fullfile(Session.spath, Session.name));
Spk.map = Map;

% SELECT specific units
if isempty(units),
    cind = true(numel(Res),1);
else
    cind = find(ismember(Clu,units));
end            
Res = Res(cind);
Clu = Clu(cind);

switch mode
% REMOVE burst tail spikes
  case 'deburst'
    nRes = [];
    nClu = [];
    thresh = round(0.008*Session.sampleRate);
    if thresh==0,thresh =1;end
    for u = unique(Clu)'
        try                            
            tRes   = Res(Clu==u);
            bursts = SplitIntoBursts(tRes,thresh);
            nRes   = [nRes; tRes(bursts)];
            nClu   = [nClu; u.*ones(numel(bursts),1)];
        end
    end
    Res = nRes;
    Clu = nClu;
  otherwise
    % NOTHING
end

% RESAMPLE to target sampleRate
Res = Res/Session.sampleRate*Spk.sampleRate;
if Spk.sampleRate~=1
    Res = ceil(Res);
end

[Res, ind] = SelectPeriods(Res,ceil(Session.sync([1,end])*Spk.sampleRate+1),'d',1,1);
Clu = Clu(ind);

% SELECT specific states
if ~isempty(states);
    if ischar(states),
        [Res,sind] = SelectPeriods(Res,[Session.stc{states,Spk.sampleRate}.data],'d',1,0);
    else
        sst = states.copy;
        sst.resample(Spk.sampleRate);
        [Res,sind] = SelectPeriods(Res,sst.data,'d',1,0);                   
    end
    Clu = Clu(sind);
end


Spk.clu = Clu;
Spk.res = Res;

Spk.update_hash(sampleRate,mode);
