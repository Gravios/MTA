function gen_pfs_figs(Trial,varargin)
[trialName,pfs_set] = DefaultArgs(varargin,{'all','basic'});

if ~isa(Trial,'MTATrial'),
    Trial = MTATrial(Trial,{},trialName);
end

switch pfs_set,
  case 'basic'
    bhvs = {'walk','rear','head','nrhp'};
  case 'basic.theta'
    bhvs = {'walk.theta','rear.theta','head.theta','nrhp.theta'};
  otherwise
    error('Unrecognized Place Field Set');
end

Pfs = {};

Trial = Trial.load_Pfs();

pf_search = MTAPlaceField([]);
pf_search.trackingMarker = 'head_front';
pf_search.smooth = 0.03;

for i = 1:numel(bhvs),
    pf_search.stateLabel = bhvs{i};
    Pfs{i} = Trial.getPfs(pf_search);
end

pfview(Trial,Pfs,trialName,'report')
