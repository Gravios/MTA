Trial = MTATrial('jg05-20120317');
lfp = Trial.lfp.copy;
lfp.create(Trial,71);
phs = lfp.phase;
%rphs = GetSegs(phs.data,round(Trial.stc{'r',lfp.sampleRate}(:,1)-3*phs.sampleRate),round(6*phs.sampleRate),0);
%figure,imagesc(rphs(:,:)'),axis  xy,

uphs = phs.copy;
Trial.load('xyz');
Trial.load('ang');
Trial.ufr.create(Trial,Trial.xyz,[],1:95,.2);
uphs.resample(Trial.ufr)

%ruphs = GetSegs(uphs.data,round(Trial.stc{'r'}(:,1)-3*uphs.sampleRate),round(6*uphs.sampleRate),0);

ow = 1;
p_t = MTAApfs(Trial,1:95,'theta',ow,'tt101',[30,30],[1,1],'xy');
p_r = MTAApfs(Trial,1:95,'rear&theta',ow,'tr101',[30,30],[1,1],'xy');
p_w = MTAApfs(Trial,1:95,'walk&theta',ow,'tw101',[30,30],[1,1],'xy');

p_t = MTAApfs(Trial,1:95,'theta',ow);
p_r = MTAApfs(Trial,1:95,'rear&theta',ow);
p_w = MTAApfs(Trial,1:95,'walk&theta',ow);

[accg,tibn] = autoccg(Trial,[],'theta',320,25,'hz','deburst');
[wccg,tibn] = autoccg(Trial,[],'walk&theta',320,25,'hz','deburst');
[rccg,tibn] = autoccg(Trial,[],'rear&theta',320,25,'hz','deburst');

units = 1:95;
f = figure;
for unit = units;
subplot(321),p_t.plot(unit);  title(p_t.parameters.states)
subplot(322),bar(tibn,accg(:,unit)),axis tight,title(num2str(unit))
subplot(323),p_r.plot(unit);  title(p_r.parameters.states)
subplot(324),bar(tibn,rccg(:,unit)),axis tight
subplot(325),p_w.plot(unit);  title(p_w.parameters.states)
subplot(326),bar(tibn,wccg(:,unit)),axis tight
saveas(f,['C:\Users\justi_000\Documents\figures\pfs_state_autoccg\' Trial.filebase '.psac-' num2str(unit) '.png']);
end



Trial.spk.create(Trial,Trial.xyz.sampleRate,'theta');

spk_r = Trial.spk.copy;
spk_r.create(Trial,Trial.xyz.sampleRate,'rear&theta',[],'deburst');
spk_w = Trial.spk.copy;
spk_w.create(Trial,Trial.xyz.sampleRate,'walk&theta',[],'deburst');
spk_t = Trial.spk.copy;
spk_t.create(Trial,Trial.xyz.sampleRate,'theta',[],'deburst');

figure,


f = figure;

units = 1:95;
for unit = units;
try
subplot(331),p_t.plot(unit);  title(p_t.parameters.states)
subplot(332),circ_plot(phs(spk_t(unit)),'hist',[],30,true,true),title(num2str(unit))
subplot(333),bar(tibn,accg(:,unit)),axis tight

subplot(334),p_r.plot(unit);  title(p_r.parameters.states)
subplot(335),circ_plot(phs(spk_r(unit)),'hist',[],30,true,true)
subplot(336),bar(tibn,rccg(:,unit)),axis tight

subplot(337),p_w.plot(unit);  title(p_w.parameters.states)
subplot(338),circ_plot(phs(spk_w(unit)),'hist',[],30,true,true)
subplot(339),bar(tibn,wccg(:,unit)),axis tight

saveas(f,['C:\Users\justi_000\Documents\figures\pfs_state_spkphs\' Trial.filebase '.pssp-' num2str(unit) '.png']);
end
end



[accg,tbins] = autoccg(Trial);

rccg = gen_bhv_ccg(Trial,'rear');
wccg = gen_bhv_ccg(Trial,'walk');

u = 91;
%figure,
for u = 1:95,
    try
subplot2(6,3,[1,2],1);p_t.plot(u);
subplot2(6,3,[1,2],2);p_r.plot(u);
subplot2(6,3,[1,2],3);p_w.plot(u);
subplot2(6,3,[3,4],1);circ_plot(uphs(spk_t(u)),'hist',[],30,false,true);
subplot2(6,3,[3,4],2);circ_plot(uphs(spk_r(u)),'hist',[],30,false,true);
subplot2(6,3,[3,4],3);circ_plot(uphs(spk_w(u)),'hist',[],30,false,true);
subplot2(6,3,[5,6],1);bar(tbins,accg(:,u));axis tight
subplot2(6,3,5,2);rccg.plot(u,1);axis tight
subplot2(6,3,6,2);rccg.plot(u,2);axis tight
subplot2(6,3,5,3);wccg.plot(u,1);axis tight
subplot2(6,3,6,3);wccg.plot(u,2);axis tight

    catch err
        continue
    end
waitforbuttonpress    
end
    