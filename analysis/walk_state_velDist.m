

s = MTASession('jg05-20120317',[],1);
s.ang.create(s);
s.ang.save;


Trial = MTATrial('jg05-20120317');
Trial.xyz.load(Trial);
Trial.filter('xyz',gausswin(11)./sum(gausswin(11)));
v = MTADxyz([],[],Filter0(gausswin(21)./sum(gausswin(21)),[zeros(1,9);Trial.vel]),Trial.xyz.sampleRate);
Trial.ang.load(Trial);


cons = [1,2,3,4,5;2,3,4,5,7];
mars = [1,3,4,7];

figure
ind = 1;
bper = Trial.stc{'t',Trial.xyz.sampleRate}.data;
for i = cons
for j = mars
subplotfit(ind,size(cons,2)*numel(mars));
hist2([log10(abs(v(bper,j))+.01),Trial.ang(bper,i(1),i(2),2)],100,100);
ind = ind+1;
end
end

cons = [1,2,3,4,5;2,3,4,5,7];
ssts = ['wtr'];

figure
ind = 1;
for j = ssts
bper = Trial.stc{j,Trial.xyz.sampleRate}.data;
for i = cons
subplotfit(ind,size(cons,2)*numel(ssts));
hist2([Trial.ang(bper,1,2,2),Trial.ang(bper,i(1),i(2),2)],100,100);
ind = ind+1;
end
end

figure,
bper = Trial.stc{'w',Trial.xyz.sampleRate}.data;
hist2([Trial.xyz(bper,5,3),Trial.ang(bper,1,3,2)],100,100);

lv = log10(abs(v(:,3)));
xaind = Trial.xyz(:,1,1)~=0&~isnan(Trial.ang(:,1,2,2))&~isinf(lv);
[sts,hmm,dcd] = gausshmm([Trial.xyz(xaind,5,3),Trial.ang(xaind,1,3,2),lv(xaind)],5);

fsts = zeros(Trial.xyz.size(1),1);
fsts(xaind) = sts;
figure,
plot(fsts),
Lines(Trial.stc{'w'}(:,1),[],'k');
Lines(Trial.stc{'w'}(:,2),[],'r');
ylim([-1,6])


%% This is the good stuff
Trial.ang.load(Trial);
lv = v.copy;
lv.data = log10(abs(lv(:,3)));
bper = Trial.stc{'w',Trial.xyz.sampleRate}.copy;
bper.cast('TimeSeries');
xaind = Trial.xyz(:,1,1)~=0&~isnan(Trial.ang(:,1,2,2))&~isinf(lv(:))&bper(:);
fet = [Trial.xyz(xaind,5,3),Trial.ang(xaind,1,3,2),lv(xaind)];
[bsts,bhmm,bdcd] = gausshmm(fet,2);


fasts = zeros(Trial.xyz.size(1),1);
fasts(xaind) = bsts;

figure,
plot(fasts),
Lines(Trial.stc{'w'}(:,1),[],'m');
Lines(Trial.stc{'w'}(:,2),[],'m');
ylim([-1,3])



Trial.stc.addState(Trial.spath,...
                   Trial.filebase,...
                   fasts==2,...
                   Trial.xyz.sampleRate,...
                   Trial.xyz.sync.copy,...
                   Trial.xyz.origin,...
                   'hwalk','g');


Trial.stc.addState(Trial.spath,...
                   Trial.filebase,...
                   fasts==1,...
                   Trial.xyz.sampleRate,...
                   Trial.xyz.sync.copy,...
                   Trial.xyz.origin,...
                   'lwalk','l');




Trial.xyz.data= [];
Trial.xyz.load(Trial),
Trial.resync(Trial.xyz),
figure,
hold on
plot([1:Session.xyz.size(1)]+Session.xyz.origin,Session.xyz(:,7,3))
plot([1:Trial.xyz.size(1)]+Trial.xyz.origin,Trial.xyz(:,7,3),'r')





sper = tper.copy;
sper.cast('TimeSeries');
Trial.resync(tper);
figure,
plot(sper.data),ylim([-1,3])
Lines(tper(:,1),[],'g');
Lines(tper(:,2),[],'g');

