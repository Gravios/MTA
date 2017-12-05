function [svar,states,stateSize,velMean,velStd] = pfs_overdispertion(Trial,varargin)

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct(stcMode...
);
[] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



%% Get units which are of good enough quality
units = select_units(Trial,25,'pyr');
Trial.load('nq');
units = units(Trial.nq.SNR(units)>1);
numClu = numel(units);

display = false;
smpleRate = 20;


xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,1,'low');
vxy = xyz.vel(Trial.trackingMarker,[1,2]);
vxy.resample(sampleRate);
xyz.resample(sampleRate);
xyz.data = sq(xyz(:,Trial.trackingMarker,[1,2]));

% LOAD unit firing rate
ufr = Trial.ufr.copy;
ufr.create(Trial,xyz,[],units);
ufr.resample(xyz);


sxyz = xyz.copy;
sufr = ufr.copy;

Tag = '';
switch mode
    case 'std'
        states = {'theta','rear','walk','hswalk','lswalk'};
% $$$     case 'rnd'
% $$$         Trial.stc.states{end+1} = rndState(Trial,'t',0.5);
% $$$         states = {'theta','rear','walk','hswalk','lswalk','x'};
% $$$         MTAApfs(Trial,units,states{end},true,'binDims',[30,30],'SmoothingWeights',[1.8,1.8]);
% $$$     case 'rnd1'
% $$$       Trial.stc.states{end+1} = rndState(Trial,'t',0.5);
% $$$       a=clock;
% $$$       Tag = num2str(a(end)+a(end-1)*100);
% $$$         states = {'x'};
% $$$         MTAApfs(Trial,units,states{end},true,Tag,'binDims',[30,30],'SmoothingWeights',[1.8,1.8]);
end
expr = {};
ufrwd ={};
pfc = {};
ethresh =5;
for s = 1:numel(states),
    sxyz.data = xyz(stc{states{s}},:);
    sufr.data = ufr(stc{states{s}},:);
    velMean = vxy(stc{states{s}});
    velMean = nanmean(log10(velMean(nniz(velMean))));
    velStd = vxy(stc{states{s}});
    velStd = nanstd(log10(velStd(nniz(velStd))));

    % LOAD placefields
    pfc{s} = MTAApfs(Trial,units,states{s},false,'binDims',[20,20],'SmoothingWeights',[1.8,1.8]);

    twpmr = ones(sxyz.size(1),numel(units));
    [~,indx] = min(abs(repmat(pfc{s}.adata.bins{1}',sxyz.size(1),1)...
                       -repmat(sxyz(:,1),1,numel(pfc{s}.adata.bins{1}))),[],2);
    [~,indy] = min(abs(repmat(pfc{s}.adata.bins{2}',sxyz.size(1),1)...
                       -repmat(sxyz(:,2),1,numel(pfc{s}.adata.bins{2}))),[],2);
    indrm = sub2ind(pfc{s}.adata.binSizes',indx,indy);


    for unit = units,
        rateMap = pfc{s}.plot(unit);
        twpmr(:,unit==units) = rateMap(indrm)/sxyz.sampleRate;
    end

    chunksz = xyz.sampleRate*ethresh;
    trim = mod(size( twpmr,1),chunksz);
    wpmr = reshape(twpmr(1:end-trim,:),chunksz,[],numel(units));
    wufr = reshape(sufr.data(1:end-trim,:),chunksz,[],numel(units));


    expr{s} = sq(sum(wufr))/sxyz.sampleRate;

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
    stateSize(s) = sum(diff(stc{states{s}.data,1,2));
end
%Trial.stc.states(end) = [];

% $$$ 
% $$$ % how does the rate varience change with state sampleSize 
% $$$ function nstate = rndState(Trial,baseStateKey,binWidth)
% $$$ RandStream.setGlobalStream(RandStream('mt19937ar','seed',sum(100*clock)));
% $$$ xyz = Trial.load('xyz');
% $$$ binSize = round(binWidth*xyz.sampleRate);
% $$$ tper = stc{baseStateKey};
% $$$ tper.cast('TimeSeries');
% $$$ tper.resample(xyz);
% $$$ tind = find(tper.data);
% $$$ newStateSize = sum(diff(stc{'w'}.data,1,2));
% $$$ nSamps = round(newStateSize/binSize);
% $$$ sampGrp = tind(randi(numel(tind),nSamps,1));
% $$$ ntind = unique(reshape(repmat(sampGrp,[1,binSize]) + repmat(1:binSize,[numel(sampGrp),1]),[],1));
% $$$ 
% $$$ nstate = tper.copy;
% $$$ nstate.data = zeros(size(nstate.data));
% $$$ nstate.data(ntind) = 1;
% $$$ nstate.cast('TimePeriods');
% $$$ nstate.label = ['rnd' tper.label];
% $$$ nstate.key = 'x';



















