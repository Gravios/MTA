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
Res = Res/Session.sampleRate*Spk.sampleRate;
if Spk.sampleRate~=1
    Res = ceil(Res);
end

[Res, ind] = SelectPeriods(Res,ceil(Session.sync([1,end])*Spk.sampleRate+1),'d',1,1);
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
    if ischar(states),
        [Res,sind] = SelectPeriods(Res,[Session.stc{states,Spk.sampleRate}.data],'d',1,0);
    else
        sst = states.copy;
        sst.resample(Spk.sampleRate);
        [Res,sind] = SelectPeriods(Res,sst.data,'d',1,0);                   
    end
    Clu = Clu(sind);
end

%% Extra Modifications
switch mode
  case 'deburst'
    nRes = [];
    nClu = [];
    thresh = round(10/1000*sampleRate);
    if thresh==0,thresh =1;end
    for u = unique(Clu)'
        try                            
            tRes = Res(Clu==u);
            bresind = SplitIntoBursts(tRes,thresh);
            nRes = [nRes; tRes(bresind)];
            nClu = [nClu; u.*ones(numel(bresind),1)];
        end
    end
    Res = nRes;
    Clu = nClu;

end


Spk.clu = Clu;
Spk.res = Res;

