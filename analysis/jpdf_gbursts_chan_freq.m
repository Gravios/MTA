
Trial = MTATrial('jg05-20120310');
stc_mode = 'auto_wbhr';
channels = 65:96;
Trial.stc.updateMode(stc_mode);Trial.stc.load;
ds = load(fullfile(Trial.spath,[Trial.name '.SelectBursts3.lfpinterp.all.1-96.mat']));
numsts = numel(Trial.stc.states);

figure(32342)
skeys = Trial.stc.list_state_attrib('key');
%skeys = {'t','r','g','l'};
skeys = cell2mat(skeys);
numsts = numel(skeys);
fbins = linspace(29,190,33);
for s = skeys,
    sper = Trial.stc{s,Trial.lfp.sampleRate};
    [~,bind] = SelectPeriods(round(ds.BurstTime*Trial.lfp.sampleRate-sper.origin),sper.data,'d',1,0);
    binds = false([numel(ds.BurstTime),1]);binds(bind)=true;
    bchanc=[];
    for i = channels,
        bchanc(:,end+1) = histc(ds.BurstFreq(ismember(ds.BurstChan,i)&binds),fbins);
    end
    subplot2(1,numsts,1,find(s==skeys));imagesc(fbins,channels,bchanc'/sum(bchanc(:)));title(Trial.stc{s}.label);xlabel('Burst Frequency');
    if s==1,ylabel('Channel'),end
end
suptitle('JPDF of detected bursts between channel and frequency given a behavioral or neural state');

lfp = Trial.lfp.copy;
lfp.load(Trial,80);% LM Phase? or pyrmidal phase?
lfp = lfp.phase;

fbins = linspace(29,190,15);
fedgs = [fbins(1:end-1);fbins(2:end)];
pedgs = linspace(-pi,pi,13);
skeys = Trial.stc.list_state_attrib('key');

j=4;
chans = [1:j:32];
figure(12325)
skeys = {'t','r','g','l'};
numsts = numel(skeys);
for s = 1:numsts,
sper = Trial.stc{skeys{s},lfp.sampleRate};
[~,bind] = SelectPeriods(round(ds.BurstTime*lfp.sampleRate-sper.origin),sper.data,'d',1,0);
binds = false([numel(ds.BurstTime),1]);binds(bind)=true;
pchanc = zeros([length(pedgs),length(fedgs)]);
for c = chans
for f = fedgs,
    btms = round(ds.BurstTime(ismember(ds.BurstChan,[64:(64+j-1)]+c)&binds&ds.BurstFreq>f(1)&ds.BurstFreq<=f(2))*lfp.sampleRate-sper.origin);
    btms(btms<0|btms>lfp.size(1))=[];
    if numel(btms)<3,continue,end
    if ~isempty(btms), pchanc(:,f(1)==fedgs(1,:)) = histc(lfp(btms),pedgs);end
end
subplot2(numel(chans),numsts,find(c==chans),s);imagesc(pedgs,fbins,(pchanc./repmat(sum(pchanc,2),1,size(pchanc,2)))');
%title(Trial.stc.states{s}.label);xlabel('Burst Frequency');
end
end


