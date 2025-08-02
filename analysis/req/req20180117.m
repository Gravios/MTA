% req20180117 ----------------------------------------------------
%  Status: active
%  Type: Analysis
%  Final_Forms: NA
%  Description: Quick Gamma Burst Analysis
%  Bugs: NA


Trial = MTATrial.validate('jg05-20120310.cof.all');
Trial = MTATrial.validate('jg05-20120311.cof.all');

stc = Trial.load('stc','msnn_ppsvd_raux');

chans = 65:96;
Trial.lfp.filename = [Trial.name,'.lfp'];
lfp = Trial.load('lfp',chans);
phz = lfp.phase([6,12]);

ds = load(['/storage/gravio/data/project/general/',Trial.name,'/',...
           Trial.name,'.DetectGammaBursts3.lfpinterp.65-96.mat']);
ds.BurstTime = ds.BurstTime-Trial.sync(1);
rmBursts = ds.BurstTime<=0|ds.BurstTime>=lfp.sync.sync.data(end);
dsFieldNames = fieldnames(ds);
dsFieldNames(cellfun(@isempty,regexp(dsFieldNames,'^Burst'))) = [];
for f = dsFieldNames',  ds.(f{1})(rmBursts)=[];  end
ds.BurstTime = round(ds.BurstTime*1250)+1;


burstFreqs = unique(ds.BurstFreq);


figure,
cs = 65:3:96;
sts = 'rxpms';
nsts = numel(sts);
for s = 1:nsts,
    for k = 1:numel(cs),
        c = cs(k):cs(k)+2;
        ind = ismember(ds.BurstChan,c)&WithinRanges(ds.BurstTime,[stc{['t&',sts(s)],lfp.sampleRate}.data]);
        subplot2(nsts,numel(cs),s,k);
        hist2([ds.BurstFreq(ind),phz(ds.BurstTime(ind),6)],burstFreqs(2:end)-diff(burstFreqs)/2,linspace(-pi,pi,20));
        title(num2str([c]))
        caxis([0,20]);
        if k ==1, ylabel(stc{sts(s)}.label);end
    end
end


