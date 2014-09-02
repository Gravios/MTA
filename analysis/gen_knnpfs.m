function gen_knnpfs(Trial,varargin)
%function gen_knnpfs(SessionName,varargin)
%
%  Ex: gen_knnpfs(Trial,true,{'rear&theta','walk&theta'})
%
%    DefaultArgs: 
%      overwrite, false
%      states,    {'theta','rear&theta','walk&theta','hwalk&theta','lwalk&theta'}
%      stc_mode   'auto_wbhr'
%      unit_type  'pyr'
%
%(varargin,{false,{'theta','rear&theta','walk&theta',,'auto_wbhr','pyr'});

[overwrite,states,stc_mode,unit_type] = DefaultArgs(varargin,{false,{'theta','rear&theta','walk&theta','hang&theta','lang&theta'},'auto_wbhr','pyr'});

%Trial = MTATrial(SessionName,TrialName);
xyz = Trial.xyz.copy; xyz.load(Trial);
Trial.load('nq');
Trial.stc.updateMode(stc_mode);Trial.stc.load;

units = select_units(Trial,18,unit_type);

if ~iscell(states), states = {states}; end
numsts = numel(states);

for i = 1:numsts,
    MTAAknnpfs(Trial,units,states{i},true,'numIter',1,...
                  'ufrShufBlockSize',0,'binDims',[10,10],...
                  'distThreshold',125,'nNearestNeighbors',110);
    %MTAAknnpfs(Trial,units,states{i},overwrite,'numIter',1000, ...
    %            'ufrShufBlockSize',0.5,'binDims',[30,30],'distThreshold',70);
end

% $$$ pfr =  MTAAknnpfs(Trial,units,states{i},true,'numIter',1000, ...
% $$$                 'ufrShufBlockSize',0.5,'binDims',[30,30],'distThreshold',70);
% $$$ 
% $$$ states = {'theta','rear&theta','walk&theta','hang&theta','lang&theta'};
% $$$ pfs ={};
% $$$ for i = 1:numel(states),
% $$$ pfs{i} = MTAAknnpfs(Trial,units,states{i},false,'numIter',1,'ufrShufBlockSize',0,'binDims',[10,10],'distThreshold',125,'nNearestNeighbors',110);
% $$$ end
% $$$ 
% $$$ figure,
% $$$ for u = units,
% $$$     for i = 1:numel(states),
% $$$         subplot(4,1,i)
% $$$         rmap = pfs{i}.plot(u);
% $$$     end
% $$$ 
% $$$     for i = 1:numel(states),
% $$$         subplot(4,1,i)
% $$$         pfs{i}.plot(u);
% $$$     end
% $$$     waitforbuttonpress
% $$$ end
% $$$ 
% $$$ 
