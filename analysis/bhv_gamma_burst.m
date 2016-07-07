
'GammaBurstPower', 'lfpinterp',  'RUN', '65,96'

ds = load(fullfile(Trial.spath,[Trial.name,'.DetectGammaBursts3_ALL.lfpinterp.65-96.mat']));

[~,tind] = SelectPeriods(ds.BurstTime,Trial.sync.data,'d',1,0);

for field = fieldnames(ds)'
    field = field{1};
    if ~isstruct(ds.(field))
        ds.(field) = ds.(field)(tind);
    end
end



ds.BurstTime = ds.BurstTime-Trial.sync(1);

s = 'w';
channels = 65:96;
fbins = linspace([ds.Params.FreqRange+[10,-10],33]);

sper = Trial.stc{s,1};
[~,bind] = SelectPeriods(ds.BurstTime,sper.data,'d',1,0);
binds = false([numel(ds.BurstTime),1]);
binds(bind)=true;
bchanc=[];
for i = channels,
    bchanc(:,end+1) = histc(ds.BurstFreq(ismember(ds.BurstChan,i)&binds),fbins);
end



figure,imagesc(fbins,channels,bsxfun(@rdivide,bchanc,sum(bchanc))')
figure,imagesc(fbins,channels,bchanc')

xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,2.4,'low');
vxy = xyz.vel(1,[1,2]);
vxy.data(vxy.data<1e-3) = 1e-4;
vxy.data = log10(vxy.data);
vxy.data = vxy(
vbins = linspace(-3,2,20);
vbinds = histc(

bchanc(:,end+1) = histc(ds.BurstFreq(ismember(ds.BurstChan,i)& binds),fbins);

bchanc=[];
for i = channels,
    bchanc(:,end+1) = histc(ds.BurstFreq(ismember(ds.BurstChan,i)),fbins);
end
figure,imagesc(fbins,channels,bsxfun(@rdivide,bchanc,sum(bchanc))')