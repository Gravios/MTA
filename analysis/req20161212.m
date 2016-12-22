'er02-20110906'
'er02-20110907'
'er02-20110908'
'er02-20110913'


Trial = MTATrial.validate(['er02-20110913','.cof.all']);

try
    stc = Trial.load('stc','NN0317');
    for s = 1:numel(stc.states),
        sts = stc.states{s};
        sts.resample(Trial.lfp.sampleRate);
        sts.data = sts.data+round(Trial.sync.data(1).*Trial.lfp.sampleRate);
        Save2sts(sts.data,fullfile(Trial.spath,[Trial.name,'.sts.',sts.label]));
    end
end

try
    stc = Trial.load('stc','NN0317R');
    for s = 1:numel(stc.states),
        sts = stc.states{s};
        sts.resample(Trial.lfp.sampleRate);
        sts.data = sts.data+round(Trial.sync.data(1).*Trial.lfp.sampleRate);
        Save2sts(sts.data,fullfile(Trial.spath,[Trial.name,'.sts.',sts.label,'_r']));
    end
end
