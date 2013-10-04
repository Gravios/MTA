
%% Load Session if Session is not already a MTASession
if ~isa(Session,'MTASession'),
    s = MTASession('jg05-20120317');
    Trial = MTATrial(s,...
                    {{'CluRes',s.xyzSampleRate}
                      'Pfs',   {'walk','rear'}},...
                    'all');
    clear('s');
end

%% Filter xyz 
Trial = Trial.filter();

%% Kernal with klen number of bins
kern = ones(klen,1);

%% number of xbins and ybins in placefields

xbins = pfw.xbin);
ybins = fliplr(pfw.ybin);

%% Get the mean position of each bin in the xy coordinates
%% Resample xyz to new bin intervals
myxyz = sq(Trial.xyz(:,7,[1,2]));
t =         permute(reshape(          myxyz(1:size(myxyz,1)-mod(size(myxyz,1),klen),:),klen,[],2),[4,1,2,3]);
for shift = 1:klen/overlap-1,
t = cat(1,t,permute(reshape(circshift(myxyz(1:size(myxyz,1)-mod(size(myxyz,1),klen),:),-overlap*shift),klen,[],2),[4,1,2,3]));
end
myxyz = t;
myxyz = reshape(sq(sum(repmat(permute(repmat(permute(repmat(kern./sum(kern),1,size(myxyz,3)),[5,4,1,2,3]),size(myxyz,4),1),[2,3,4,1]),klen/overlap,1).*myxyz,2)),[],2);


%% Calculate the new sampling rate for indexing purposes 
newSampleRate = 1/((size(Trial.xyz,1)-mod(size(Trial.xyz,1),klen))/Trial.xyzSampleRate/length(myxyz));


%% Selecte the xy periods from the down-sampled 
%% UNUSED AT THE MOMENT
% $$$ myrx = SelectPeriods(myxyz,round((Trial.Bhv.getState('rear').state+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),'c',1,1);
% $$$ mywx = SelectPeriods(myxyz,round((Trial.Bhv.getState('walk').state+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),'c',1,1);
% $$$ myrv = sqrt(sum(diff(myrx).^2,2));
% $$$ mywv = sqrt(sum(diff(mywx).^2,2));
