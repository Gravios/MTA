function frameSet_series(Trial,SetName,skip)

Trial= MTATrial('jg05-20120317');
SetName = 'frameset_852_1438';
SetName = 'frameset_243_447';



load(fullfile(Trial.spath, [Trial.filebase '.' SetName '.mat']));
imPath =fullfile(Trial.spath, [Trial.filebase '.' SetName '/']);
if ~exist(imPath,'dir'),
    mkdir(imPath);
end

indices = [];
images = {};
for i = 1:length(record),
    indices(i) = record{i}.index;
    record{i} = rmfield(record{i},'index');
    images{i} = frame2im(record{i});
end

%% browse frames

hfig = figure;
nsp = 5;
skip = 10;
fcount = 9;
graythresh = 230;graygrad = round(10:60/fcount:60);gcount = 0;
inds = round(linspace(1,numel(images),numel(images)/skip));
setname = 'grayed_rear_sequence';
mov = [];

qflip = 0;
qmov = 1;

swrp = exp(-linspace(-1,1,nsp).^2);
swrpn = swrp./sum(swrp);

for i = 1:nsp
    sp(i) = subplot(1,nsp,i);
end   
for i = 1:nsp
    pos = get(sp(i),'outerposition');
    pos(1) = sum(swrpn(1:i))-swrpn(i);
    pos(2) = (1-swrp(i))./2;
    pos(3) = swrpn(i);
    pos(4) = swrp(i);
    set(sp(i),'outerposition',pos);
end


h=[];hi=[];api=[];
index = ceil(nsp/2);
for i = 1:nsp,
    hi(i) = image(images{inds(index-floor(nsp/2)+i)},'Parent',sp(i));
end

while index ~=-1;
    for i = 1:nsp,       
        set(hi(i),'CData',images{inds(index-floor(nsp/2)+i)},'Parent',sp(i));
    end
    B = waitforbuttonpress;
    whatkey = get(hfig,'CurrentCharacter');
    if ~B,continue,end
    switch double(whatkey)
        case double('n')
            index = index+1;    
        case double('p')
            index = index-1;    
        case double('q')
            index = -1;
        case double('s')
            if ishandle(h),delete(h);end
            [h,api]=imRectRot('rotate',0,'hParent',sp(round(nsp/2)));
            uistack(h,'top');
            cpos = round(feval(api.getPos));
        case double('c')
            cim = get(hi(ceil(nsp/2)),'CData');
            croppedIm = cim(cpos(2):(cpos(4)+cpos(2)),cpos(1):(cpos(3)+cpos(1)),:);
            if qflip,croppedIm = fliplrim(croppedIm);end
            imwrite(croppedIm,fullfile(imPath,[ setname '-cropped_frame_' num2str(indices(inds(index-floor(nsp/2)+ceil(nsp/2)))) '.png']),'PNG','BitDepth',16,'Transparency',[1,1,1]);
            if qmov,mov = cat(4,mov,croppedIm);end
            index = index+1;    
        case double('g')
            gcount = gcount+1;
            cim = get(hi(ceil(nsp/2)),'CData');
            croppedIm = cim(cpos(2):(cpos(4)+cpos(2)),cpos(1):(cpos(3)+cpos(1)),:);
            grim = repmat(mean(croppedIm,3),[1,1,3]);
            grim = grim./max(grim(:));
            contrastRange = [graythresh-graygrad(gcount),255];
            x = round(diff(contrastRange).*grim)+contrastRange(1);
            x(x(:)==min(x(:)))=1;
            x = uint8(x);
            if qflip,x = fliplrim(x);end
            imwrite(x,fullfile(imPath,[ setname '-cropped_frame_' num2str(indices(inds(index-floor(nsp/2)+ceil(nsp/2)))) '.png']),'PNG','BitDepth',16,'Transparency',[1,1,1]);%,'Alpha', ones(size(x,1),size(x,2)).*alphagrad(gcount));
            if qmov,mov = cat(4,mov,x);end
            index = index+1;    
    end
    if (index+floor(nsp/2)+1)>=numel(inds),index = 2;end
    if index<ceil(nsp/2)&&index~=-1,index = numel(inds)-ceil(nsp/2);end
end

fshift = 100;
F = filmStrip( mov(:,:,:,1:1:9), size(mov,2)-fshift, 0, 0 ,1);
figure,image(F)


