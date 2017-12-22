function frameSet_viewer(Trial,SetName,skip)

Trial= MTATrial('jg05-20120317');
SetName = 'frameset_852_1438';
skip = 20;


load(fullfile(Trial.spath, [Trial.filebase '.' SetName '.mat']));
imPath =fullfile(Trial.spath, [Trial.filebase '.' SetName '/']);
if ~exist(imPath,'dir'),
    mkdir(imPath);
end


index = [];
images = {};
for i = 1:length(record),
index(i) = record{i}.index;
record{i} = rmfield(record{i},'index');
images{i} = frame2im(record{i});
end



nid = numel(images);
figure
for i = linspace(1,nid,round(nid/skip)),
    image(images{i});
end
