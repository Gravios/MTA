
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
uind = nq.Refrac<5e-4;
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
uind = nq.Refrac<5e-4;
figure, hold on
plot(mpfs.patchMFR(s1,uind,1),...
     mpfs.patchMFR(s2,uind,1),...
     '.')
% $$$ plot(mpfs.patchPFR(s1,uind,1),...
% $$$      mpfs.patchPFR(s2,uind,1),...
% $$$      '.')
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


% TEST VARS ---------------------------------------------------------------------------
Trial = MTATrial.validate('jg05-20120309.cof.all');
unit = 27;
unit = 85;
state = 'pause&theta'
state = 'loc&theta'
%--------------------------------------------------------------------------------------




% HELPER Functions --------------------------------------------------------------------

% GET the trial index within the sesList struct array
gti = @(sesList,Trial) find(arrayfun(@(s,t) strcmp([s.sessionName,'.',s.mazeName,'.',s.trialName],t{1}.filebase),...
                                     sesList,...
                                     repmat({Trial},size(sesList))));
% GET the state index within the states cell array
gsi = @(states,state) find(cellfun(@strcmp,states,repmat({state},size(states))),1);
%--------------------------------------------------------------------------------------




% MAIN --------------------------------------------------------------------------------

tind = gti(sesList,Trial);
suind = find(map(:,4)==tind);
units = map(suind,1);

sper = Trial.stc{state};
sper.cast('TimeSeries');
sper.resample(rhm);


[~,mpind] = max(mpfs.patchMFR(sind,suind(units==unit),:));
%mpfs.patchRateInd(sind,suind(units==unit),mpind,:,:)
% $$$ figure,
% $$$ plot(sq(mpfs.patchRateInd(sind,suind(units==unit),mpind,1,:)),...
% $$$      sq(mpfs.patchRateInd(sind,suind(units==unit),mpind,2,:)),'.');


% Trial specific stuff
sind = gsi(states,state);                            % get state index
Trial.load('stc',pfstats{tind}.stcMode);             % load (neural network/hand) labeled states
[rhm,fs,ts] = fet_rhm(Trial,[],'mtchglong',true);    % compute rhythmic head motion power
ufr = Trial.ufr.copy;                                % load unit firing rates
ufr.create(Trial,rhm,'theta-groom-sit',units,0.8);
xyz = Trial.load('xyz');                             % load xyz coordinates
xyz.resample(rhm);                                   % resample xyz to rhm sampleRate
%ufr.resample(rhm);                                   % resample ufr to rhm sampleRate
pft = pfs_2d_theta(Trial);                           % load placefields during theta periods

% xyz position mapped onto place field bins' indicies
[~,indx] = min(abs( repmat(pft.adata.bins{1}',xyz.size(1),1)...
                    -repmat(xyz(:,Trial.trackingMarker,1),1,numel(pft.adata.bins{1}))),...
               [],2);
[~,indy] = min(abs( repmat(pft.adata.bins{2}',xyz.size(1),1)...
                    -repmat(xyz(:,Trial.trackingMarker,2),1,numel(pft.adata.bins{2}))),...
               [],2);
rpow = MTADxyz('data',log10(median(rhm(:,6<fs&fs<13),2)),'sampleRate',rhm.sampleRate);



% select time indicies where xyz position is within a place field's best patch
patchInd = and(ismember(indx,sq(mpfs.patchRateInd(sind,suind(units==unit),mpind,1,:))),...
               ismember(indy,sq(mpfs.patchRateInd(sind,suind(units==unit),mpind,2,:))));

ind = patchInd&sper.data==1;

figure,plot(rpow(ind),log10(ufr(ind,units==unit)+abs(randn([sum(ind),1]))./10),'.')

figure,plot(rpow(ind),ufr(ind,units==unit),'.')

% rhm distribution
eds = linspace(-8.5,-2,100);
figure,bar(eds,histc(rpow(nniz(rpow)),eds),'histc');


% NEXT bin mean ufr binned by rpow
