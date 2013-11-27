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
comv = cat(1,comv(1,:),comv);


%% body center of mass
%% no filter
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


%% center of mass log10 speed
%% distribution of rbb rigid body
%% during walking vs everything else
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% body center of mass
%% filter - gausswin(121)
combvf = Filter0(gausswin(121)./sum(gausswin(121)),comv(:,1))*12;
figure,
windowlen = 40000;
windowoffset = 1;
plot(log10(clip(combvf(windowoffset:windowoffset+windowlen,1),0.0001,50)))
xlabel(['time in frame @' num2str(Trial.xyzSampleRate)]);
ylabel(['com body log10 speed cm/s']);
Lines([],0,'k');
wper = Trial.Bhv.getState('walk').state;
Lines(wper(wper(:,1)<windowlen&wper(:,1)>windowoffset,1),[],'g');
Lines(wper(wper(:,1)<windowlen&wper(:,1)>windowoffset,2),[],'m');
rper = Trial.Bhv.getState('rear').state;
Lines(rper(rper(:,1)<windowlen&rper(:,1)>windowoffset,1),[],'k');
Lines(rper(rper(:,1)<windowlen&rper(:,1)>windowoffset,2),[],'r');
reportfig([], 'ex20130417', 0,...
['center of mass log10 speed of ' num2str(cell2mat(rbb.ml)) ' rigid body with walking ' 'onset(g) offset(m) and rear onset(k) offset(r)']);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% center of mass filtered (gausswin(121))
%% log10 speed distribution of rbb rigid body
%% during walking vs everything else
figure
bhvl = 'walk'
bhvdist = SelectPeriods(combvf(:,1),Trial.Bhv.getState(bhvl).state,'c',1,1);
notbhvper = SubstractRanges([1,size(combvf,1)],Trial.Bhv.getState('walk').state);
notbhvdist = SelectPeriods(combvf(:,1),notbhvper,'c',1,1);
edges = [-1.4:7/1000:2];
nnwd = histc(log10(abs(notbhvdist)),edges);
nwd  = histc(log10(abs(bhvdist)),edges);
figure,
subplot(211)
bar(edges,nwd,'histc')
ylabel('count')
xlabel(['com body log10 speed cm/s']);
yl = ylim;
subplot(212)
bar(edges,nnwd,'histc')
ylabel('count')
xlabel(['com body log10 speed cm/s']);
ylim(yl.*[1,2])
reportfig([], 'ex20130417', 0,...
['center of mass log10 speed distribution of ' num2str(cell2mat(rbb.ml)) ...
 ' rigid body with walking ' 'onset(g) offset(m) vs everything else'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% center of mass filtered (gausswin(121))
%% log10 speed distribution of rbb rigid body
%% during rearing vs everything else
figure
bhvl = Trial.Bhv.list_States;
yl = zeros(length(bhvl),2);
for i = 1:length(bhvl),
bhvdist = SelectPeriods(combvf(:,1),Trial.Bhv.getState(bhvl{i}).state,'c',1,1);
notbhvper = SubstractRanges([1,size(combvf,1)],Trial.Bhv.getState(bhvl{i}).state);
notbhvdist = SelectPeriods(combvf(:,1),notbhvper,'c',1,1);
edges = [-1.4:7/1000:2];
nnwd = histc(log10(abs(notbhvdist)),edges);
nwd  = histc(log10(abs(bhvdist)),edges);
subplot(length(bhvl)+1,1,i)
bar(edges,nwd,'histc')
title(bhvl{i});
ylabel('count')
xlabel(['com body log10 speed cm/s']);
yl(i,:) = ylim;
end
subplot(length(bhvl)+1,1,i+1)
bar(edges,nnwd,'histc')
ylabel('count')
xlabel(['com body log10 speed cm/s']);
ylm = max(yl(:));
ylim([0,ylm])

reportfig([], 'ex20130417', 0,...
['center of mass log10 speed distribution of ' num2str(cell2mat(rbb.ml)) ...
 ' rigid body with rearing ' 'onset(g) offset(m) vs everything else'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


