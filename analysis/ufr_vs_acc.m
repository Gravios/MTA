%% UFR vs ACC

win = 241;
comb = Filter0(gausswin(win)./sum(gausswin(win)),Trial.com(Trial.Model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'})));
comh = Filter0(gausswin(win)./sum(gausswin(win)),Trial.com(Trial.Model.rb({'head_back','head_left','head_front','head_right'})));

combv = sqrt(sum(diff(comb(:,[1,2])).^2,2));
comhv = sqrt(sum(diff(comh(:,[1,2])).^2,2));

ufrsegs = {};
chvsegs = {};
for i = 1:length(group_restrains),
ufrsegs{i} = GetSegs(Trial.ufr(:,16),round((group_restrains{i}-1)./Trial.xyzSampleRate.*Trial.lfpSampleRate)-2*Trial.lfpSampleRate,4*Trial.lfpSampleRate,0);
chvsegs{i} = GetSegs(abs(diff(comhv)),round(group_restrains{i}-2*Trial.xyzSampleRate),round(4*Trial.xyzSampleRate),0);
end
ufrslen = size(ufrsegs{1},1);
chvslen = size(chvsegs{1},1);

for i = 1:length(group_restrains),
ufrsegs{i} = ufrsegs{i}(round(Trial.lfpSampleRate/Trial.xyzSampleRate/2:Trial.lfpSampleRate/Trial.xyzSampleRate:ufrslen),:);
end


figure,plot(ufrsegs{6},chvsegs{6},'.')
figure,plot(sum(ufrsegs{6},2),sum(chvsegs{6},2),'.')

for i = 1:length(group_restrains),
figure,plot(mean(ufrsegs{i},2),mean(chvsegs{i},2),'.')
end
