load ~/data/analysis/jg05-20120317/jg05-20120317.cof.all.bayes2.mat

s = MTASession('jg05-20120317');

Trial = MTATrial(s,{{'CluRes',s.xyzSampleRate}},'all');

Trial = Trial.filter();

klen = 16;
kern = ones(klen,1);
overlap = 4;
nnn = 81;
nxbins = 50;
nybins = 50;
xbins = linspace(Trial.Maze.boundaries(1,1),Trial.Maze.boundaries(1,2),nxbins)';
ybins = linspace(Trial.Maze.boundaries(2,2),Trial.Maze.boundaries(2,1),nxbins)';


myxyz = sq(Trial.xyz(:,7,[1,2]));
t =         permute(reshape(          myxyz(1:size(myxyz,1)-mod(size(myxyz,1),klen),:),klen,[],2),[4,1,2,3]);
for shift = 1:klen/overlap-1,
t = cat(1,t,permute(reshape(circshift(myxyz(1:size(myxyz,1)-mod(size(myxyz,1),klen),:),-overlap*shift),klen,[],2),[4,1,2,3]));
end
myxyz = t;
myxyz = reshape(sq(sum(repmat(permute(repmat(permute(repmat(kern./sum(kern),1,size(myxyz,3)),[5,4,1,2,3]),size(myxyz,4),1),[2,3,4,1]),klen/overlap,1).*myxyz,2)),[],2);

newSampleRate = 1/((size(Trial.xyz,1)-mod(size(Trial.xyz,1),klen))/Trial.xyzSampleRate/length(myxyz));

myrx = SelectPeriods(myxyz,round((Trial.Bhv.getState('rear').state+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),'c',1,1);
mywx = SelectPeriods(myxyz,round((Trial.Bhv.getState('walk').state+0.5*Trial.xyzSampleRate/newSampleRate)./Trial.xyzSampleRate.*newSampleRate),'c',1,1);


% $$$ unit = 1;
% $$$ ufrs = {};
% $$$ %%iu = [12,14,24];
% $$$ while unit ~=-1


for unit = 1:76,
myres = Trial.res(Trial.clu==unit);
myufr{unit} = zeros(size(Trial.xyz,1),1);
myufr{unit}(1:myres(end)) = accumarray(myres,ones(size(myres),1));
t =         permute(reshape(myufr{unit}(1:size(myufr{unit},1)-mod(size(myufr{unit},1),klen)),klen,[]),[3,1,2]);
for shift = 1:klen/overlap-1
t = cat(1,t,permute(reshape(circshift(myufr{unit}(1:size(myufr{unit},1)-mod(size(myufr{unit},1),klen)),-overlap*shift),klen,[]),[3,1,2]));
end
myufr{unit} = t; 
myufr{unit} = reshape(sq(sum(repmat(permute(repmat(kern,1,size(myufr{unit},3)),[3,1,2]),size(myufr{unit},1),1).*myufr{unit},2))/(klen/Trial.xyzSampleRate),[],1);
try,ufrs{unit} = unique(myufr{unit});end
end

ufr = cell2mat(myufr);

tic
dcdim = zeros(size(ufr,1),76);
for unit = 1:76,
for i = 1:size(ufr,1),
dcdim(i,unit) = find(ufrs{unit}==ufr(i,unit));
end
end
toc

% $$$ pmap = zeros(50,50,76,size(ufr,1));
% $$$ for unit = 9:76,
% $$$ for i = 1:size(ufr,1),
% $$$     %if dcdim(i,unit)~=1,
% $$$ pmap(:,:,unit,i) = pfknnmrw{unit}(:,:,dcdim(i,unit));
% $$$ %end
% $$$ end
% $$$ end


Smooth = 0.03;
msize = [50,50];
r1 = (-msize(1):msize(1))/msize(1);
r2 = (-msize(2):msize(2))/msize(2);
Smoother1 = exp(-r1.^2/Smooth^2/2);
Smoother1 = Smoother1/sum(Smoother1);
Smoother2 = exp(-r2.^2/Smooth^2/2);
Smoother2 = Smoother2/sum(Smoother2);




pxyz = zeros(size(myxyz));
figure
for i = 1:size(myxyz,1),
    clf
    pmap =[];
    pmapr =[];
    for unit = gu,
        if dcdim(i,unit)<0,continue,end
        onan = pfknnmrw{unit}(:,:,dcdim(i,unit));
        lnan = pfknnmrw{unit}(:,:,dcdim(i,unit));
        lnan(isnan(lnan))=0;
        %np = conv2(Smoother1,Smoother2,lnan,'same');
        np = lnan;
        np(isnan(onan)) = 1;
        if sum(np(:))==0,continue,end
        pmap(:,:,end+1) =np;

% $$$     onanr = pfknnmrr{unit}(:,:,dcdim(i,unit));
% $$$     lnanr = pfknnmrr{unit}(:,:,dcdim(i,unit));
% $$$     lnanr(isnan(lnanr))=0;
% $$$     %np = conv2(Smoother1,Smoother2,lnan,'same');
% $$$     npr = lnanr;
% $$$     npr(isnan(onanr)) = 1;
% $$$     if sum(npr(:))==0,continue,end
% $$$     pmapr(:,:,end+1) =npr;
end
if sum(pmap(:)==0)==numel(pmap),pmap = ones(50,50,length(gu));,end
pmap = pmap(:,:,2:end);
%pmapr = pmapr(:,:,2:end);
pmap = cat(3,pmap,pmapr);
predmap = -1./log10(prod(pmap,3));
predmap(isinf(predmap))=nan;
%predmap = prod(pmap,3);
%predmap = sum(pmap,3);
bw =  predmap>prctile(predmap(:),99);
bw = imfill(bw, 'holes');
L = bwlabel(bw);
s = regionprops(L, 'PixelIdxList', 'PixelList');
if size(s,1)~=0
idx = cat(1,s.PixelIdxList);
sum_region = sum(predmap(idx));
pixelList = cat(1,s.PixelList);
x = pixelList(:, 1);
y = pixelList(:, 2);
xbar = round(sum(x .* double(predmap(idx))) / sum_region);
ybar = round(sum(y .* double(predmap(idx))) / sum_region);
else
xbar = 25;
ybar = 25;
end
imagescnan({xbins,ybins,predmap},[],[],0,[0,0,0]);
hold on,
scatter(xbins(xbar),ybins(ybar),20,[1,0,0],'fill')
scatter(myxyz(i,1),myxyz(i,2),20,[1,1,1],'fill')
pause(0.01)
if ~isnan(xbar)&~isnan(ybar)
pxyz(i,:) = [xbins(xbar),ybins(ybar)];
end
xlim([-500,500])
ylim([-500,500])
end



pfw = MTAPlaceField(Trial,[],'walk');
gu = [];
for index =1:76,
pfw.plot(index),
waitforbuttonpress
whatkey = get(gcf,'CurrentCharacter');
switch double(whatkey)
  case double('n')
    index = index+1;
  case double('m')
    gu(end+1) = index;
  case double('i')
    index = input('Enter index #: ');
  case double('p')
    index=index-1;
  case double('q')
    index = -1;
end
end
