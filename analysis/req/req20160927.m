% req20160927

clear
TrialList = {'Ed01-20140707.cof.all','Ed03-20140624.cof.all', ...
             'Ed03-20140625.cof.all','Ed05-20140529.ont.all'};

hfet = {};
hfetCorrected = {};
xyz = {};
vxy = {};
for t = TrialList
    RefTrial = MTATrial.validate(TrialList{2});    
    Trial = MTATrial.validate(t{1});
    xyz {end+1} = Trial.load('xyz');
    xyz{end}.filter('ButFilter',3,2.4,'low');
    vxy {end+1} = xyz{end}.vel([1,7],[1,2]);
    hfet{end+1} = fet_head_pitch(Trial);
    hfetCorrected{end+1} = hfet{end}.copy;
    hfetCorrected{end}.map_to_reference_session(Trial,RefTrial);
end


alpha = 0.5;
ct = 'rgcm';
eds = linspace(-pi/2,pi/2,200);
figure;
subplot(211);
hold('on');
for t = 1:numel(TrialList),
    %ind = nniz(hfet{t});    
    ind = nniz(hfet{t})&vxy{t}(:,1)>3;
    hs = bar(eds,histc(hfet{t}(ind),eds),'histc');
    hs.FaceColor = ct(t);
    hs.EdgeColor = ct(t);
    hs.FaceAlpha = alpha;
    hs.EdgeAlpha = alpha;
end
subplot(212);
hold('on');
for t = 1:numel(TrialList),
    %ind = nniz(hfetCorrected{t});
    ind = nniz(hfet{t})&vxy{t}(:,1)>3;
    hs = bar(eds,histc(hfetCorrected{t}(ind),eds),'histc');
    hs.FaceColor = ct(t);
    hs.EdgeColor = ct(t);
    hs.FaceAlpha = alpha;
    hs.EdgeAlpha = alpha;
end




mid = kmeans(mzd(nniz(mzd)),2);
alpha = 0.5;
ct = 'rgcm';
eds = linspace(-pi/2,pi/2,200);
figure;
subplot(211);
hold('on');
ind = nniz(mzd);
val = mzd(ind);    

for t = 1:2
    hs = bar(eds,histc(val(mid==t),eds),'histc');
    median(val(mid==t))
    hs.FaceColor = ct(t);
    hs.EdgeColor = ct(t);
    hs.FaceAlpha = alpha;
    hs.EdgeAlpha = alpha;
end
