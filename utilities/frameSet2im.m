function frameSet2im(Session,SetName)

SetName = 'frameset_16768_17797';
load([Session.spath.analysis Session.filebase '.' SetName '.mat'])
imPath =[Session.spath.analysis Session.filebase '.' SetName '/'];
if ~exist(imPath,'dir'),
    mkdir(imPath);
end

index = [];
images = {};
for i = 1:length(record),
index(i) = record{i}.index;
record{i} = rmfield(record{i},'index');
images{i} = frame2im(record{i});
%imwrite(images{i},[imPath 'frame_' num2str(index(i)) '.png'],'PNG','BitDepth',16,'Transparency',[0,0,0])
end


cimages = {}
for i = 1:length(record),
cimages{i} = images{i}(1:520,200:750,:);
end


figure,
for i = 1:length(record),
imwrite(cimages{i},[imPath 'frame_crop' num2str(index(i)) '.png'],'PNG','BitDepth',16)
end
