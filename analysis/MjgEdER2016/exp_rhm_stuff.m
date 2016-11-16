

% LOAD DATA ---------------------------------------------------------------------------
% Summary
%   N := number of units
%   S := number of states
%
%   map    - Nx4 matrix with unit identification
%   nq     - Nx1 struct with array fields regarding unit characteristics   
%   mpfs   - SxN struct with array fields regarding place field characteristics
%   states - 1xS cell array of strings with behavior state labels
%            corresponding to all other varialbles with an S dim

load_pfs_unit_stats;

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
