

Trial = MTATrial('jg05-20120310');
fet = detect_walk(Trial);
figure,plot(fet)
Lines(Trial.stc{'w'}(:),[],'k');
Lines(Trial.stc{'r'}(:),[],'w');
Lines([],0,'k');

fet = MTADxyz('data',fet,'sampleRate',Trial.xyz.sampleRate);






sfet = sign(fet.data);
afet = abs(fet.data);
afet(afet<1)=1;
lfet = log10(afet).*sfet;

fet.data = lfet;

figure,hist(fet(Trial.stc{'r'}),100)

Trial.xyz.filter(gausswin(61)./sum(gausswin(61)));
vel = Trial.vel;
vel = clip(log10(vel),-2,2);
vel =  MTADxyz('data',[zeros([1,size(vel,2)]);vel],'sampleRate',Trial.xyz.sampleRate);


m=1;
hist2([vel(:,m),fet(:)],linspace(-2,2,64),linspace(-2,2,64));
caxis([0,600])



figure
hist2([vel(:,7),clip(log10(Trial.xyz(:,7,3)),0,3)],linspace(-2,2,64),linspace(1.6,2.6,64))
caxis([0,1600])
hold on
vh = [vel(Trial.stc{'r'},7),clip(log10(Trial.xyz(Trial.stc{'r'},7,3)),0,3)];
eeh = error_ellipse(cov(vh),mean(vh),'conf',.95,'style','-w');
vh = [vel(Trial.stc{'w'},7),clip(log10(Trial.xyz(Trial.stc{'w'},7,3)),0,3)];
eeh = error_ellipse(cov(vh),mean(vh),'conf',.95,'style','-g');
vh = [vel(Trial.stc{'l'},7),clip(log10(Trial.xyz(Trial.stc{'l'},7,3)),0,3)];
eeh = error_ellipse(cov(vh),mean(vh),'conf',.95,'style','-c');
vh = [vel(Trial.stc{'g'},7),clip(log10(Trial.xyz(Trial.stc{'g'},7,3)),0,3)];
eeh = error_ellipse(cov(vh),mean(vh),'conf',.95,'style','-m');


figure
ind = ':';
vm = 1;
hm = 1;
hist2([vel(:,vm),clip(log10(Trial.xyz(:,hm,3)),0,3)],linspace(-2,2,64),linspace(1,2,64))
caxis([0,1600])
hold on
states = 'rwgl';
colors = 'wgmc';
for i = 1:numel(states),
ind = Trial.stc{states(i)};
vh = [vel(ind,vm),clip(log10(Trial.xyz(ind,hm,3)),0,3)];
eeh = error_ellipse(cov(vh),mean(vh),'conf',.95,'style',colors(i));
end


figure
ind = ':';
vm = 1;
hm = 1;
hist2([vel(:,vm),clip(log10(Trial.xyz(:,hm,3)./10),0,3)],linspace(-2,2,64),linspace(0.2,1,64))
ticks_lin2log
caxis([0,1600])
hold on
states = 'rwgl';
colors = 'wgmc';
for i = 1:numel(states),
ind = Trial.stc{states(i)};
vh = [vel(ind,vm),clip(log10(Trial.xyz(ind,hm,3)./10),0,3)];
eeh = error_ellipse(cov(vh),mean(vh),'conf',.95,'style',colors(i));
end


fetset = fetset_generic(Trial);

ind = 1700;
figure,imagesc(unity(fetset(ind:ind+800,1:8))')
Lines(Trial.stc{'w'}(:)-ind,[],'k');
Lines(Trial.stc{'r'}(:)-ind,[],'b');



figure
hist2([clip(log10(Trial.xyz(:,1,3)),0,3),clip(log10(Trial.xyz(:,7,3)),0,3)],linspace(1,2,64),linspace(1.6,2.6,64))
%caxis([0,1600])
hold on
vh = [clip(log10(Trial.xyz(Trial.stc{'r'},1,3)),0,3),clip(log10(Trial.xyz(Trial.stc{'r'},7,3)),0,3)];
eeh = error_ellipse(cov(vh),mean(vh),'conf',.95,'style','-w');
vh = [clip(log10(Trial.xyz(Trial.stc{'w'},1,3)),0,3),clip(log10(Trial.xyz(Trial.stc{'w'},7,3)),0,3)];
eeh = error_ellipse(cov(vh),mean(vh),'conf',.95,'style','-g');
vh = [clip(log10(Trial.xyz(Trial.stc{'l'},1,3)),0,3),clip(log10(Trial.xyz(Trial.stc{'l'},7,3)),0,3)];
eeh = error_ellipse(cov(vh),mean(vh),'conf',.95,'style','-c');
vh = [clip(log10(Trial.xyz(Trial.stc{'g'},1,3)),0,3),clip(log10(Trial.xyz(Trial.stc{'g'},7,3)),0,3)];
eeh = error_ellipse(cov(vh),mean(vh),'conf',.95,'style','-m');


figure
hist2([clip(log10(Trial.xyz(:,1,3)),0,3),clip(log10(Trial.xyz(:,7,3)),0,3)],linspace(1,2,64),linspace(1.6,2.6,64))
%caxis([0,1600])
hold on
vh = [clip(log10(Trial.xyz(Trial.stc{'r'},1,3)),0,3),clip(log10(Trial.xyz(Trial.stc{'r'},7,3)),0,3)];
eeh = error_ellipse(cov(vh),mean(vh),'conf',.95,'style','-w');
vh = [clip(log10(Trial.xyz(Trial.stc{'w'},1,3)),0,3),clip(log10(Trial.xyz(Trial.stc{'w'},7,3)),0,3)];
eeh = error_ellipse(cov(vh),mean(vh),'conf',.95,'style','-g');
vh = [clip(log10(Trial.xyz(Trial.stc{'l'},1,3)),0,3),clip(log10(Trial.xyz(Trial.stc{'l'},7,3)),0,3)];
eeh = error_ellipse(cov(vh),mean(vh),'conf',.95,'style','-c');
vh = [clip(log10(Trial.xyz(Trial.stc{'g'},1,3)),0,3),clip(log10(Trial.xyz(Trial.stc{'g'},7,3)),0,3)];
eeh = error_ellipse(cov(vh),mean(vh),'conf',.95,'style','-m');





hist2([vel(:,7),clip(log10(Trial.xyz(:,'head_front',3)),1,3)],100,100)

hist2([clip(log10(v(:,'head_front')),-2,3),clip(log10(Trial.xyz(1:end-1,'head_front',3)),1.6,3)],100,100)

hist2([clip(log10(v(:,'head_front')),-2,3),clip(log10(Trial.xyz(1:end-1,'head_front',3)),1.6,3).*Trial.ang(1:end-1,3,