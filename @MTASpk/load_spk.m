function Spk = load_spk(Spk,Session,varargin)
% function Spk = create(Spk,Session,varargin)
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
%     loadField: cellarry, Auxillary options for LoadCluResSpk
%
% See also MTAStateCollection for more information on selecting
% states
%
%% Load and resample Res             
defargs = struct('sampleRate',               1,   ...
                 'states',                   [],...
                 'units',                    [],...
                 'mode',                     '',...
                 'loadField',                 {{}} ...
                 );
[sampleRate,states,units,mode,loadField] = DefaultArgs(varargin,defargs,'--struct');
Spk.sampleRate = sampleRate;
[clu,res,spk,map] = LoadCluResSpk(fullfile(Session.spath, Session.name),loadField{:});

Spk.map = map;
res = res/Session.sampleRate*Spk.sampleRate;
if Spk.sampleRate~=1
    res = ceil(res);
end

[res, ind] = SelectPeriods(res,ceil(Session.sync([1,end])*Spk.sampleRate+1),'d',1,1);
clu = clu(ind);
spk = spk(ind,:,:);

%% Select specific units
if isempty(units),
    cind = true(numel(res),1);
else
    cind = find(ismember(clu,units));
end            

res = res(cind);
clu = clu(cind);
spk = spk(cind,:,:);

%% Select specific states
if ~isempty(states);
    if ischar(states),
        [res,sind] = SelectPeriods(res,[Session.stc{states,Spk.sampleRate}.data],'d',1,0);
    else
        sst = states.copy;
        sst.resample(Spk.sampleRate);
        [res,sind] = SelectPeriods(res,sst.data,'d',1,0);                   
    end
    clu = clu(sind);
    spk = spk(sind,:,:);
end

%% Extra Modifications
switch mode
  case 'deburst'
    nRes = [];
    nClu = [];
    nSpk = [];
    thresh = round(10/1000*sampleRate);
    if thresh==0,thresh =1;end
    for u = unique(clu)'
        try                            
            tRes = res(clu==u);
            tSpk = spk(clu==u,:,:);
            bresind = SplitIntoBursts(tRes,thresh);
            nRes = [nRes; tRes(bresind)];                        
            nClu = [nClu; u.*ones(numel(bresind),1)];
            nSpk = cat(1,nSpk,tSpk(bresind,:,:));
        end
    end
    res = nRes;
    clu = nClu;
    spk = nSpk;

end


Spk.clu = clu;
Spk.res = res;
Spk.spk = spk;            
   