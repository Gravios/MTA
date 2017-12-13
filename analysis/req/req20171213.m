


sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);
pitchReferenceTrial = 'Ed05-20140529.ont.all';

Trials  = af(@(t)  MTATrial.validate(t),   sessionList);
units   = cf(@(t)  select_placefields(t),  sessionList);

numTrials = numel(Trials);

pfd  = cell([1,3]);  
[pfd{:}]= cf(@(t)  MjgER2016_drzfields(t), Trials);
highRateInds = -0.5 < pfdPitch{1}.adata.bins{1} & pfdPitch{1}.adata.bins{1} < 0.5;

rateMapInd = nan([numTrials,max(cellfun('length',units)),numel(pfd)]);
rateMapVal = nan([numTrials,max(cellfun('length',units)),numel(pfd)]);
rateMaxRng = nan([numTrials,max(cellfun('length',units)),numel(pfd)],2);
for t = 1:numTrials,
    for u = 1:numel(units{t}),
        for p = 1:numel(pfd),
            rateMap = plot(pfd{p}{t},units{t}(u),'isCircular',false);
            rateMap = mean(rateMap(highRateInds,:),'omitnan');
            [rateMaxVal(t,u,p),rateMaxInd(t,u,p)] = max(rateMap);
            rateMaxRng(t,u,p,:) = ThreshCross(rateMap,rateMaxVal(t,u,p)/2,1);
        end
    end
end

rateMaxVal = sq(reshape(rateMaxVal,[],1,numel(pfd)));
rateMaxInd = sq(reshape(rateMaxInd,[],1,numel(pfd)));
rateMaxRng = sq(reshape(rateMaxRng,[],1,numel(pfd),2));

figure();

binsPitch  = pfd{1}{1}.adata.bins{2};
binsHeight = pfd{2}{1}.adata.bins{2};
binsRHM    = pfd{2}{1}.adata.bins{2};

rateMaxPitch  = binsPitch(rateMaxInd(:,1));
rateMaxHeight = binsHight(rateMaxInd(:,2));
rateMaxRHM    = binsRHM(rateMaxInd(:,3));


        