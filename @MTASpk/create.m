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
if isempty(units)
    cind = true(numel(Res),1);
elseif isnumeric(units)
    cind = find(ismember(Clu,units));
elseif ischar(units)
    cind = find(ismember(Clu,get_unit_set(Spk,Session,units)));
end            
Res = Res(cind);
Clu = Clu(cind);

switch mode
% REMOVE burst tail spikes
  case 'deburst'
    nRes = [];
    nClu = [];
    thresh = round(0.008*Session.sampleRate);
    for u = unique(Clu)'
        try                            
            tRes   = Res(Clu==u);
            burstResInd = SplitIntoBursts(tRes,thresh);
            nRes   = [nRes; tRes(burstResInd)];
            nClu   = [nClu; u.*ones(numel(burstResInd),1)];
        end
    end
    Res = nRes;
    Clu = nClu;
  case 'blge2'
    nRes = [];
    nClu = [];
    thresh = round(0.008*Session.sampleRate);
    blThresh = 2;    
    for u = unique(Clu)'
        try
            tRes   = Res(Clu==u);
            [burstResInd,burstLength] = SplitIntoBursts(tRes,thresh);
            nRes   = [nRes; tRes(burstResInd(burstLength>=blThresh))];
            nClu   = [nClu; u.*ones([sum(burstLength>=blThresh),1])];
        end        
    end
    Res = nRes;
    Clu = nClu;    
    
  case 'blge3'
    nRes = [];
    nClu = [];
    thresh = round(0.008*Session.sampleRate);
    blThresh = 3;
    for u = unique(Clu)'
        try
            tRes   = Res(Clu==u);
            [burstResInd,burstLength] = SplitIntoBursts(tRes,thresh);
            nRes   = [nRes; tRes(burstResInd(burstLength>=blThresh))];
            nClu   = [nClu; u.*ones([sum(burstLength>=blThresh),1])];
        end        
    end
    Res = nRes;
    Clu = nClu;    
  case 'first_spike_theta'
    % get first spike of each theta cycle
    
    
  otherwise % default
    % NOTHING
end

% RESAMPLE to target sampleRate
Res = Res/Session.sampleRate*Spk.sampleRate;
if Spk.sampleRate~=1
    Res = ceil(Res);
end


if ~isempty(Spk.per)
    newRes = [];
    newClu = [];
    for cluInd = unique(Clu)'
        sper = copy(Spk.per);
        sper.data = sper.data(Spk.perInd(cluInd,:),:)./sper.sampleRate.*Spk.sampleRate;
        if isempty(sper.data)
            continue
        end
        [pRes,psind] = SelectPeriods(Res,sper.data,'d',1,0);
        newRes = cat(1, newRes, pRes);
        newClu = cat(1, newClu, cluInd*ones(size(pRes)));
    end
    [Res,sind] = sort(newRes);
    Clu        = newClu(sind);
end


[Res, ind] = SelectPeriods(Res,ceil(Session.sync([1,end])*Spk.sampleRate+1*double(Spk.sampleRate~=1)),'d',1,1);
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
