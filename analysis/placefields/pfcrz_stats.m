function pfcrz_stats(batch_trial_list)

batch_trial_list = ['/data/homes/gravio/data/analysis/' ...
                    'batch_trial_list'];
trialName = 'all';

zs = [];
rs = [];
ms = [];

fid = fopen(batch_trial_list);
SessionName = fgets(fid);
while ischar(SessionName)
    SessionName = SessionName(1:end-1);
    Trial = MTATrial(SessionName,[],trialName);
    pfcrz = MTAPlaceField(Trial,[],'theta',0,trialName,'pfcrz');
    pfw = MTAPlaceField(Trial,[],'walk');
    numClu = size(pfw.cluMap,1);
    for unit = 1:numClu
        if isempty(pfw.maxRateMax{unit}),continue,end
        if isempty(pfcrz.maxRateMax{unit}),continue,end
        if pfw.maxRate{unit}(pfw.maxRateMax{unit})<5,continue,end
        if pfcrz.maxRate{unit}(pfcrz.maxRateMax{unit})<2,continue,end
        ms(end+1) = nanmedian(pfw.rateMap{unit}(:));
        if ms(end)>1.4,continue,end
        rs(end+1) = pfcrz.maxRatePos{unit}(pfcrz.maxRateMax{unit},1);
        zs(end+1) = pfcrz.maxRatePos{unit}(pfcrz.maxRateMax{unit},2);
    end
    SessionName = fgets(fid);
end
fclose(fid)
