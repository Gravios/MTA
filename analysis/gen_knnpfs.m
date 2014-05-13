function gen_knnpfs(SessionName,varargin)
%function gen_knnpfs(SessionName,varargin)
%
%  Ex: gen_knnpfs('jg05-20120312',1,{'rear&theta','walk&theta'})
%
%    DefaultArgs: 
%      TrialName, 'all'
%      overwrite, false
%      states,    {'theta','rear&theta','walk&theta','hwalk&theta','lwalk&theta'}
%      stc_mode   'auto_wbhr'
%      unit_type  'pyr'
%
%(varargin,{'all',false,{'theta','rear&theta','walk&theta',,'auto_wbhr','pyr'});

[TrialName,overwrite,states,stc_mode,unit_type] = DefaultArgs(varargin,{'all',false,{'theta','rear&theta','walk&theta','hwalk&theta','lwalk&theta'},'auto_wbhr','pyr'});

Trial = MTATrial(SessionName,TrialName);
Trial.xyz.load(Trial);
Trial.load('nq');
Trial.stc.updateMode(stc_mode);Trial.stc.load;

units = select_units(Trial,18,unit_type);

if ~iscell(states), states = {states}; end
numsts = numel(states);

for i = 1:numsts,
    MTAAknnpfs(Trial,units,states{i},overwrite,'numIter',1000, ...
                'ufrShufBlockSize',0.5,'binDims',[20,20],'distThreshold',70);
end
