function [svar,states,stateSize] = pfs_overdispertion(Trial,mode)


if ~strcmp(Trial.stc.mode,'auto_wbhr'),Trial.stc.updateMode('auto_wbhr');Trial.stc.load;end

%% Get units which are of good enough quality
units = select_units(Trial,18,'pyr');
numClu = numel(units);

display = false;
downSampleRate = 20;
newSampleRate = downSampleRate;

xyz = Trial.load('xyz');
myxyz = xyz.copy;
myxyz.filter(gtwin(.1,xyz.sampleRate));
myxyz.resample(newSampleRate);
myxyz.data = sq(myxyz(:,Trial.trackingMarker,[1,2]));



myufr = Trial.ufr.copy;
myufr.create(Trial,xyz,[],units);
myufr.resample(myxyz);

smyxyz = myxyz.copy;
smyufr = myufr.copy;

switch mode
    case 'std'
states = {'theta','rear','walk','hswalk','lswalk'};
    case 'rnd'
        Trial.stc.states{end+1} = rndState(Trial,'t',0.5);
        states = {'theta','rear','walk','hswalk','lswalk','x'};
        MTAApfs(Trial,units,states{end},true,'binDims',[20,20],'SmoothingWeights',[1.8,1.8]);
    case 'rnd1'
        Trial.stc.states{end+1} = rndState(Trial,'t',0.5);
        states = {'x'};
        MTAApfs(Trial,units,states{end},true,'binDims',[20,20],'SmoothingWeights',[1.8,1.8]);
end
expr = {};
ufrwd ={};
pfc = {};
ethresh =5;
for s = 1:numel(states),
smyxyz.data = myxyz(Trial.stc{states{s}},:);
smyufr.data = myufr(Trial.stc{states{s}},:);


%% Get the expected ufr for each xy 
%% Substract the expected ufr from the observed
pfc{s} = MTAApfs(Trial,units,states{s},'binDims',[20,20],'SmoothingWeights',[1.8,1.8]);
twpmr = ones(smyxyz.size(1),numel(units));
[~,indx] = min(abs(repmat(pfc{s}.adata.bins{1}',smyxyz.size(1),1)-repmat(smyxyz(:,1),1,numel(pfc{s}.adata.bins{1}))),[],2);
[~,indy] = min(abs(repmat(pfc{s}.adata.bins{2}',smyxyz.size(1),1)-repmat(smyxyz(:,2),1,numel(pfc{s}.adata.bins{2}))),[],2);
indrm = sub2ind(pfc{s}.adata.binSizes',indx,indy);

for unit = units,
    rateMap = pfc{s}.plot(unit);
    twpmr(:,unit==units) = rateMap(indrm)/smyxyz.sampleRate;
end

chunksz = myxyz.sampleRate*ethresh;
trim = mod(size( twpmr,1),chunksz);
wpmr = reshape(twpmr(1:end-trim,:),chunksz,[],numel(units));
wufr = reshape(smyufr.data(1:end-trim,:),chunksz,[],numel(units));


expr{s} = sq(sum(wufr))/smyxyz.sampleRate;

ufrwd{s} = (expr{s}-sq(sum(wpmr)))./sqrt(expr{s});

end


if display,
    hfig = figure;
    set(hfig,'position',[360,377,1070,324])
    unit = units(1);
    while unit~=-1,
        u =find(units==unit);
        for s = 1:numel(states);
            subplot2(2,numel(states),1,s),
            pfc{s}.plot(unit,[],1);
            title([states{s},': ',num2str(unit)])
            subplot2(2,numel(states),2,s),
            hist(ufrwd{s}(expr{s}(:,u)>2,u),100)
            title(['var: ',num2str(var(ufrwd{s}(nniz(ufrwd{s}(:,u))&expr{s}(:,u)>ethresh,u)))])
        end
        unit = figure_controls(hfig,unit,units);
    end
end


for s =1:numel(states),
    svar(s) = var(ufrwd{s}(nniz(ufrwd{s}(:))&expr{s}(:)>ethresh));
    stateSize(s) = sum(diff(Trial.stc{states{s}}.data,1,2));
end
Trial.stc.states(end) = [];
% how does the rate varience change with state sampleSize 
function nstate = rndState(Trial,baseStateKey,binWidth)
xyz = Trial.load('xyz');
binSize = round(binWidth*xyz.sampleRate);
tper = Trial.stc{baseStateKey};
tper.cast('TimeSeries');
tper.resample(xyz);
tind = find(tper.data);
newStateSize = sum(diff(Trial.stc{'w'}.data,1,2));
nSamps = round(newStateSize/binSize);
sampGrp = tind(randi(numel(tind),nSamps,1));
ntind = unique(reshape(repmat(sampGrp,[1,binSize]) + repmat(1:binSize,[numel(sampGrp),1]),[],1));

nstate = tper.copy;
nstate.data = zeros(size(nstate.data));
nstate.data(ntind) = 1;
nstate.cast('TimePeriods');
nstate.label = ['rnd' tper.label];
nstate.key = 'x';



















