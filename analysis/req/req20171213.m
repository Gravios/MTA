


sessionListName = 'MjgER2016';
sessionList = get_session_list(sessionListName);
pitchReferenceTrial = 'Ed05-20140529.ont.all';

Trials  = af(@(t)  MTATrial.validate(t),   sessionList);
units   = cf(@(t)  select_placefields(t),  Trials);

numTrials = numel(Trials);

pfd  = cell([1,3]);  
[pfd{:}]= cf(@(t)  MjgER2016_drzfields(t), Trials);
highRateInds = -0.5 < pfd{1}{1}.adata.bins{1} & pfd{1}{1}.adata.bins{1} < 0.5;

map        = nan([numTrials,max(cellfun('length',units)),2]);
rateMapInd = nan([numTrials,max(cellfun('length',units)),numel(pfd)]);
rateMapVal = nan([numTrials,max(cellfun('length',units)),numel(pfd)]);
rateMaxRng = nan([numTrials,max(cellfun('length',units)),numel(pfd),2]);

for t = 1:numTrials,
    for u = 1:numel(units{t}),
        for p = 1:numel(pfd),
            rateMap = plot(pfd{p}{t},units{t}(u),'isCircular',false);
            rateMap = mean(rateMap(highRateInds,:),'omitnan');
            [rateMaxVal(t,u,p),rateMaxInd(t,u,p)] = max(rateMap);
            try
                rateMaxRng(t,u,p,:) = cat(4,find([rateMap>[rateMaxVal(t,u,p)/2]]==1,1,'first'),...
                                            find([rateMap>[rateMaxVal(t,u,p)/2]]==1,1,'last'));
            end
            map(t,u,:) = cat(3,t,units{t}(u));
        end
    end
end

rateMaxVal = sq(reshape(rateMaxVal,[],1,numel(pfd)));
rateMaxInd = sq(reshape(rateMaxInd,[],1,numel(pfd)));
rateMaxRng = sq(reshape(rateMaxRng,[],1,numel(pfd),2));

binsPitch  = pfd{1}{1}.adata.bins{2};
binsHeight = pfd{2}{1}.adata.bins{2};
binsRHM    = pfd{3}{1}.adata.bins{2};

rateMaxPitch  = binsPitch(rateMaxInd(nniz(rateMaxInd,1),1));
rateMaxHeight = binsHeight(rateMaxInd(nniz(rateMaxInd,2),2));
rateMaxRHM    = binsRHM(rateMaxInd(nniz(rateMaxInd,3),3));

figure();
subplot(131);  hist(rateMaxPitch,25);
subplot(132);  hist(rateMaxHeight,20);
subplot(133);  hist(rateMaxRHM,25);

figure();
subplot(131);  plot(rateMaxPitch,diff(rateMaxRng(nniz(rateMaxRng,1,1),1,:),1,3),'.');
subplot(132);  plot(rateMaxHeight,diff(rateMaxRng(nniz(rateMaxRng,2,1),2,:),1,3),'.');
subplot(133);  plot(rateMaxRHM,diff(rateMaxRng(nniz(rateMaxRng,3,1),3,:),1,3),'.');


figure();
subplot(131);  
hold('on');
for i = 1:10,
    plot(rateMaxPitch+randn(size(rateMaxPitch))/50,rateMaxHeight+randn(size(rateMaxHeight))*10,'.b');
end

subplot(132);  
hold('on');
for i = 1:10,
    plot(rateMaxPitch+randn(size(rateMaxPitch))/50,rateMaxRHM+randn(size(rateMaxRHM))/10,'.b');
end


subplot(133);
hold('on');
for i = 1:10,
    plot(rateMaxHeight+randn(size(rateMaxHeight))*10,rateMaxRHM+randn(size(rateMaxRHM))/10,'.b');
end