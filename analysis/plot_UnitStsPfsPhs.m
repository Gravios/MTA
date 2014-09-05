%Trial = MTATrial('er06-20130614','all-cof');
Trial = MTATrial('jg05-20120317');
Trial.load('xyz');
lfp = Trial.lfp.copy;
lfp.create(Trial,61);
lfp.resample(Trial.xyz);
phs = lfp.phase;


%rphs = GetSegs(phs.data,round(Trial.stc{'r',lfp.sampleRate}(:,1)-3*phs.sampleRate),round(6*phs.sampleRate),0);
%figure,imagesc(rphs(:,:)'),axis  xy,

%uphs = phs.copy;


%Trial.load('ang');
%Trial.ufr.create(Trial,Trial.xyz,[],1:95,.2);
%uphs.resample(Trial.ufr)

%ruphs = GetSegs(uphs.data,round(Trial.stc{'r'}(:,1)-3*uphs.sampleRate),round(6*uphs.sampleRate),0);

ow = 0;
p_t = MTAApfs(Trial,[],'theta',ow,'tt101',[30,30],[1,1],'xy');
p_r = MTAApfs(Trial,[],'rear&theta',ow,'tr101',[30,30],[1,1],'xy');
p_w = MTAApfs(Trial,[],'walk&theta',ow,'tw101',[30,30],[1,1],'xy');

%p_t = MTAApfs(Trial,1:95,'theta',ow);
%p_r = MTAApfs(Trial,1:95,'rear&theta',ow);
%p_w = MTAApfs(Trial,1:95,'walk&theta',ow);

[accg,tibn] = autoccg(Trial,[],'theta',320,25,'hz','deburst');
[wccg,tibn] = autoccg(Trial,[],'walk&theta',320,25,'hz','deburst');
[rccg,tibn] = autoccg(Trial,[],'rear&theta',320,25,'hz','deburst');

spk_r = Trial.spk.copy;
spk_r.create(Trial,Trial.xyz.sampleRate,'rear&theta',[],'deburst');
spk_w = Trial.spk.copy;
spk_w.create(Trial,Trial.xyz.sampleRate,'walk&theta',[],'deburst');
spk_t = Trial.spk.copy;
spk_t.create(Trial,Trial.xyz.sampleRate,'theta',[],'deburst');


%rccg = gen_bhv_ccg(Trial,'rear');
%rccg = gen_bhv_ccg(Trial,'rear');
%wccg = gen_bhv_ccg(Trial,'walk');


%f = figure;

units = 1:size(Trial.spk.map,1);
for unit = units;
try
subplot(331) ;cla; p_t.plot(unit);  title(p_t.parameters.states)
subplot(332) ;cla; circ_plot(phs(spk_t(unit)),'hist',[],30,true,true),title(num2str(unit))
subplot(333) ;cla; bar(tibn,accg(:,unit)),axis tight
catch
subplot(331);cla
subplot(332);cla
subplot(333);cla
end

try
subplot(334);cla; p_r.plot(unit);  title(p_r.parameters.states)
subplot(335);cla; circ_plot(phs(spk_r(unit)),'hist',[],30,true,true)
subplot(336);cla; bar(tibn,rccg(:,unit)),axis tight
catch
subplot(334);
subplot(335);
subplot(336);
end

try
subplot(337);cla; p_w.plot(unit);  title(p_w.parameters.states)
subplot(338);cla; circ_plot(phs(spk_w(unit)),'hist',[],30,true,true)
subplot(339);cla; bar(tibn,wccg(:,unit)),axis tight
catch
subplot(337);
subplot(338);
subplot(339);
end

%saveas(f,['C:\Users\justi_000\Documents\figures\pfs_state_spkphs\' Trial.filebase '.pssp-' num2str(unit) '.png']);
saveas(f,['/gpfs01/sirota/bach/homes/gravio/figures/pfs_state_spkphs/' Trial.filebase '.pssp-' num2str(unit) '.png']);

end






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
    