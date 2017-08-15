


sessionList = 'hand_labeled';
Trials = af(@(Trial) MTATrial.validate(Trial), get_session_list(sessionList));
stc = cf(@(t) t.load('stc'), Trials);
%stc = cf(@(t) t.load('stc','msnnN0+hand_labeled'), Trials);

s = 1
xyz = cf(@(Trial) Trial.load('xyz'),                         Trials);
vxy = cf(@(x)     x.vel({'spine_lower','head_front'},[1,2]), xyz   );
for s = 1:numel(Trials), vxy{s}.data(vxy{s}.data<1e-3,:) = 1e-3;end
for s = 1:numel(Trials), vxy{s}.data = log10(vxy{s}.data);end
cf(@(x) x.filter('ButFilter',3,5,'low'),xyz);;
ang = cf(@(t,x) create(MTADang,t,x),Trials,xyz);

s = 1;
figure,
plot(ang{s}(:,'head_back','head_front',3))

figure();
eds  = linspace(-3,3,500);
zang = cf(@(a,x) clip(a(nniz(x),'head_back','head_left',3)...
                      -median(a(nniz(x),'head_back','head_left',3)),-3.01,3.01),ang,xyz);
for s=1:numel(zang), zang{s}(zang{s}>3|zang{s}<-3)=[];end
zang = cf(@(a) clip(a-median(a),-3.01,3.01),zang);
zang = cat(1,zang{:});
bar(eds,histc(zang(1:25:end,1),eds),'histc');
Lines(-std(zang(1:25:end,1)),[],'c');
Lines(std(zang(1:25:end,1)),[],'c');
Lines(-std(zang(1:25:end,1)*2),[],'g');
Lines(std(zang(1:25:end,1))*2,[],'g');
ylabel('Count');
xlabel('Median Centered Inter-Marker Distance (mm)')
xlim([-3,3])
title(['Inter Marker std of distance: ', num2str(round(std(zang(1:25:end,1)),4))])


