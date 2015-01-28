
Trial = MTATrial('jg05-20120310');
lfp = Trial.load('lfp',80);
phs = lfp.phase([100,120]);
spk = Trial.load('spk',lfp.sampleRate,'theta');

for i = 1:100,
    res = spk(i);
    if ~isempty(res)&&numel(res)<10000,
y = PPC(phs(res))
end
end