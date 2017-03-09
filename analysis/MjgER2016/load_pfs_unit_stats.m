
% DEFARGS -----------------------------------------------------------------------------
sesList = get_session_list('MjgEdER2016_bhv');
stcMode = 'NN0317R';
states =  {'loc','rear','pause','lloc','hloc','lpause','hpause'};
states = cellfun(@strcat,states,repmat({'&theta'},size(states)),'UniformOutput',false);
%--------------------------------------------------------------------------------------




% LOAD DATA ---------------------------------------------------------------------------
% Summary
%   N := number of units
%   S := number of states
%
%   map - Nx4 matrix with unit identification
%   nq  - struct with Nx1 array fields regarding unit characteristics   
%   mpfs- struct with SxN array fields regarding place field characteristics
%

% Load all trials
Trials = arrayfun(@MTATrial.validate,sesList,'UniformOutput',false);
% Load place field stats for all trials
pfstats = cellfun(@(Trial) batch_compute_pfstats_bs(Trial),Trials,'UniformOutput',false);
% Unwrap pfstats
pfstats = cellfun(@(x) x{1},pfstats,'UniformOutput',false);

sesIds =[1:5,9,10];

% create clumap
map = cellfun(@cat,...
              mat2cell(2*ones([1,numel(sesIds)]),1,ones([numel(sesIds),1])),...
              cellfun(@(T,S) T.spk.map(S.cluMap,:),Trials(sesIds),pfstats(sesIds),'UniformOutput',false),...
              cellfun(@times,...
                      cellfun(@(S) ones([numel(S.cluMap),1]),pfstats(sesIds),'UniformOutput',false),...
                      mat2cell(sesIds,1,ones(numel(sesIds),1)),...
                      'UniformOutput',false),...
              'UniformOutput',false);
map = vertcat(map{:});

% Load unit info
cellfun(@load,Trials(sesIds),repmat({'nq'},[1,numel(sesIds)]),'UniformOutput',false);
nq = cellfun(@(T) StructArray(T.nq,1),Trials(sesIds),'UniformOutput',false);
nq = cellfun(@(N,S) N(S.cluMap),nq,pfstats(sesIds),'UniformOutput',false);
nq = CatStruct(cat(1,nq{:}),fieldnames(nq{1}),1);

% Remodel pfstats substructure data into single struct
mpfs = cellfun(@(x) x.pfmstats,pfstats(sesIds),'UniformOutput',false);
clear('rpfstats')
fnames = fieldnames(mpfs{1})';
for t = 1:numel(mpfs),
    for f = fnames
        f = f{1};
        rpfstats{t}.(f) = sq(reshape([mpfs{t}(:,:).(f)],[size(mpfs{t}),size(mpfs{t}(1,1).(f),3),size(mpfs{t}(1,1).(f),4),size(mpfs{t}(1,1).(f),5)]));
    end
end
mpfs = CatStruct(cat(1,rpfstats{:}),fnames,2);
states = pfstats{1}.states;
stcModes = cellfun(@(x) x.stcMode,pfstats,'UniformOutput',false);
%--------------------------------------------------------------------------------------


clear('rpfstats','t','f','sesIds','pfstats')