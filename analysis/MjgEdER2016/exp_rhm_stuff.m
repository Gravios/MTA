
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
%--------------------------------------------------------------------------------------




% FIGURES -----------------------------------------------------------------------------
% Plotting some Group stats
% 
% Rear x Loc
% 
% Rear x Pause
%
% Loc  x Pause
%

s1 = 1;
s2 = 2;
uind = nq.eDist>25;
figure, hold on
plot(mpfs.patchMFR(s1,uind,1).*log10(mpfs.patchArea(s1,uind,1)),...
     mpfs.patchMFR(s2,uind,1).*log10(mpfs.patchArea(s2,uind,1)),...
     '.')
plot([0,100],[0,100],'-m')
plot([0,50],[0,100],'-r')
plot([0,100],[0,50],'-r')
plot([0,25],[0,100],'-r')
plot([0,100],[0,25],'-r')
daspect([1,1,1])


s1 = 1;
s2 = 2;
uind = nq.eDist>25;
figure, hold on
plot(mpfs.patchMFR(s1,uind,1),...
     mpfs.patchMFR(s2,uind,1),...
     '.')
plot(mpfs.patchPFR(s1,uind,1),...
     mpfs.patchPFR(s2,uind,1),...
     '.')
plot([0,40],[0,40],'-m')
plot([0,20],[0,40],'-r')
plot([0,40],[0,20],'-r')
plot([0,10],[0,40],'-r')
plot([0,40],[0,10],'-r')
daspect([1,1,1])


s1 = 1;
s2 = 2;
figure,plot(mpfs.patchMFR(s1,:,1),...
            mpfs.patchMFR(s2,:,1),...
            '.')


s1 = 1;
s2 = 2;
figure,plot(log10(mpfs.peakFR(s1,uind,1)),...
            log10(mpfs.peakFR(s2,uind,1)),...
            '.')
daspect([1,1,1])
hold on,plot(log10(linspace(0.1,40,100)),log10(linspace(0.1,40,100)),'-m')
hold on,plot(log10(linspace(0.1,20,100)),log10(linspace(0.1,40,100)),'-r')
hold on,plot(log10(linspace(0.1,40,100)),log10(linspace(0.1,20,100)),'-r')


figure,hist(log10(mpfs.peakFR(s1,uind,1)./mpfs.peakFR(s2,uind,1)),100)

%--------------------------------------------------------------------------------------



% RHM versus UFR ----------------------------------------------------------------------
% testing the rhm effects on a single session before writing
% standalone function
%
% motivation - most behaviors are permeated with sniffing.
%
% technical considerations:
%  
%   Problem:
%     Units are spatially selective and their firing rates
%     depend on location. 
%   Solution:
%     Only regress firing rate onto rhm power within idividual
%     place field patches.
%    
% Independent Vars:
%
%   rhm - rhythmic head motion average spectral power within the
%   range of 6-14 Hz.
%
%   ufr - unit firing rate with 0.8 second box car average of spikes.
%

Trial = MTATrial.validate('jg05-20120309.cof.all');
pfstats = batch_compute_pfstats_bs(Trial);

% helper function to find the patch index which contains the max firing rate
mind = @(x,s,u) find(sq(x{1}.pfmstats(s,u).patchPFR)==max(x{1}.pfmstats(s,u).patchPFR));

mind(pfstats,1,9)

for i = 1:size(tpfs,1),
    for j = 1:size(tpfs,2),
        tpfs(i,j).peakFR = permute(tpfs(i,j).peakFR,[2,3,1]);
        tpfs(i,j).rateThreshold = permute(tpfs(i,j).rateThreshold,[2,3,1]);        
        tpfs(i,j).spatialCoherence = permute(tpfs(i,j).spatialCoherence,[2,3,1]);        
    end
end


peakfr = [];
for i = 1:size(tpfs,1)

end


fnames = fieldnames(tpfs)';
clear('rpfstats')
for f = fnames
    f = f{1};
    rpfstats.(f) = sq(reshape([tpfs(:,:).(f)],[size(tpfs),size(tpfs(1,1).(f),3),size(tpfs(1,1).(f),4),size(tpfs(1,1).(f),5)]));
end

mind = @(x,s) sum(bsxfun(@times,sq(bsxfun(@eq,rpfstats.patchPFR(s,:,:),max(rpfstats.patchPFR(s,:,:),[],3))),1:2),2);

s1 = 1;
s2 = 2;
figure,plot(rpfstats.patchMFR(s1,:,1).*log10(rpfstats.patchArea(s1,:,1)),...
            rpfstats.patchMFR(s2,:,1).*log10(rpfstats.patchArea(s2,:,1)),...
            '.')






