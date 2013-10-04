Trial = MTATrial('jg05-20120309',{'ang'},'all');


win = 7;
xyz = reshape(Filter0(gausswin(win)./sum(gausswin(win)),Trial.xyz),size(Trial.xyz,1),size(Trial.xyz,2),size(Trial.xyz,3));
v = sqrt(sum(diff(xyz).^2,3));

win = 241;
rbb = Trial.Model.rb({'spine_lower','pelvis_root','spine_middle'});
rbh = Trial.Model.rb({'head_back','head_left','head_front','head_right'});
comb = Filter0(gausswin(win)./sum(gausswin(win)),Trial.com(rbb));
comh = Filter0(gausswin(win)./sum(gausswin(win)),Trial.com(rbh));

comb = sqrt(sum(diff(comb(:,[1,2])).^2,2));
comh = sqrt(sum(diff(comh(:,[1,2])).^2,2));


comv = [comb,comh];
for i =1:2,
comv(:,i) = ButFilter(comv(:,i),11,2./(Trial.xyzSampleRate/2),'low');
comv(:,i) = clip(comv(:,i),0.003,30);
end

sfet = [circ_dist(Trial.ang(:,2,3,1),Trial.ang(:,1,2,1)),...
        circ_dist(Trial.ang(:,3,4,1),Trial.ang(:,2,3,1)),...
        circ_dist(Trial.ang(:,4,5,1),Trial.ang(:,3,4,1)),...
        circ_dist(Trial.ang(:,5,7,1),Trial.ang(:,4,5,1))];

num_fet = size(sfet,2);
for i = 1:num_fet
   sfet(isnan(sfet(:,i)),i) = circ_mean(sfet(~isnan(sfet(:,i)),i));
end

win = 7;
sfet = Filter0(gausswin(win)./sum(gausswin(win)),sfet);


bsfet = ButFilter(sum(sfet,2),11,6./(Trial.xyzSampleRate/2),'low');
tsfet = Filter0(gausswin(141)./sum(gausswin(141)),bsfet);


dbsfet = diff(bsfet-tsfet); 

nfet = -Filter0(gausswin(141)./sum(gausswin(141)),abs(dbsfet.*Filter0(gausswin(141)./sum(gausswin(141)),comv(:,1)).*12))./(Trial.xyz(1:end-1,7,3)./Trial.ang(1:end-1,3,4,2)).*1000;
nfet = cat(1,nfet(1),nfet);

plot(nfet)
Lines([],0,'k');
Lines(Trial.Bhv.getState('walk').state(:,1),[],'g');
Lines(Trial.Bhv.getState('walk').state(:,2),[],'m');
Lines(Trial.Bhv.getState('rear').state(:,1),[],'k');
Lines(Trial.Bhv.getState('rear').state(:,2),[],'r');


windowlen = 40000;
windowoffset = 1;
plot(log10(clip(comv(windowoffset:windowoffset+windowlen,1),0.0001,14)))
Lines([],0,'k');

wper = Trial.Bhv.getState('walk').state;
Lines(wper(wper(:,1)<windowlen&wper(:,1)>windowoffset,1),[],'g');
Lines(wper(wper(:,1)<windowlen&wper(:,1)>windowoffset,2),[],'m');
rper = Trial.Bhv.getState('rear').state
Lines(rper(rper(:,1)<windowlen&rper(:,1)>windowoffset,1),[],'k');
Lines(rper(rper(:,1)<windowlen&rper(:,1)>windowoffset,2),[],'r');

reportfig([], 'ex20130417', 0,...
['center of mass log10 speed of ' num2str(cell2mat(rbb.ml)) ' rigid body with walking ' 'onset(g) offset(m) and rear onset(k) offset(r)']);




bhvl = 'walk'
bhvdist = SelectPeriods(comv(:,1),Trial.Bhv.getState(bhvl).state,'c',1,1);
notbhvper = SubstractRanges([1,size(comv,1)],Trial.Bhv.getState('walk').state);
notbhvdist = SelectPeriods(comv(:,1),notbhvper,'c',1,1);

edges = [-3:7/1000:1];
nnwd = histc(log10(abs(notbhvdist)),edges);
nwd  = histc(log10(abs(bhvdist)),edges);
figure,
subplot(211)
bar(edges,nwd,'histc')
yl = ylim;
subplot(212)
bar(edges,nnwd,'histc')
ylim(yl.*[1,2])
reportfig([], 'ex20130417', 0,...
['center of mass log10 speed distribution of ' num2str(cell2mat(rbb.ml)) ...
 ' rigid body with walking ' 'onset(g) offset(m) vs everything else'])



hist(nfet,1000)



bhvl = 'rear'
bhvdist = SelectPeriods(nfet,Trial.Bhv.getState(bhvl).state,'c',1,1);

notbhvper = SubstractRanges([1,size(nfet,1)],Trial.Bhv.getState('walk').state);

notbhvdist = SelectPeriods(nfet,notbhvper,'c',1,1);

edges = [-6:7/1000:1];

nnwd = histc(log10(abs(notbhvdist)),edges);
nwd  = histc(log10(abs(bhvdist)),edges);



figure,
subplot(211)
bar(edges,nwd,'histc')
subplot(212)
bar(edges,nnwd,'histc')
reportfig([], 'new_rear_walk_feature', 0,...
['distributions of ' bhvl ' and non-' bhvl ' periods for the new feature'
)
figure
plot(nfet)
Lines([],0,'k');
Lines(Trial.Bhv.getState('walk').state(:,1),[],'g');
Lines(Trial.Bhv.getState('walk').state(:,2),[],'m');
Lines(Trial.Bhv.getState('rear').state(:,1),[],'k');
Lines(Trial.Bhv.getState('rear').state(:,2),[],'r');
Lines([],10^-1.5,'k');
xlim