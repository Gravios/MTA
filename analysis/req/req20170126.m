% NAME : req20170126 
% GOAL : Examine current inter-session/subject feature coregistration

RefTrial = MTATrial.validate('jg05-20120317.cof.all');

Trial = MTATrial.validate('Ed03-20140624.cof.all');

fet = fet_mis(Trial);
rfet = fet_mis(RefTrial);

dbstop in map_to_reference_session at 126

fet.map_to_reference_session(Trial,RefTrial);


f = 1;


figure
subplot(2,2,1);
imagesc(tarMean{f}');

figure,
z = 5:10:50;
for i = 1:5,
    subplot(1,5,i);
    imagesc(tarCnt{f}(:,:,z(i))');
end

f = 1;
embeddedFeatureMapConvolution = convn(tarCnt{f},refCnt{f},'same');

mxy = [];
cxy = [];
for x = 1:25,
    for y = 1:25,
        [mxy(x,y),cxy(x,y)] = max(conv(sq(tarCnt{f}(x,y,:)),sq(refCnt{f}(x,y,:)),'same'));
    end
end
cxy = cxy(:);
mxy = mxy(:);
figure,hist2([log10(mxy(mxy~=0)),cxy(mxy~=0)],5,25);caxis([0,100])
tarDom{f}(round(mean(cxy(mxy>10))))-tarDom{f}(25)

figure,hist(tarDom{f}(cxy(mxy>10))-tarDom{f}(25))

figure,hist(cxy(nnz&mxy0))


[fsc] = max(sq(reshape(embeddedFeatureMapConvolution,[],50)));
[~,fsci] = max(fsc);
tarDom{f}(fsci)-tarDom{f}(25)

i = 25;
surf(1:50,1:50,i,cz{1}(:,:,i))
%surf(bin3(i)*ones(pfSize),rateMap(:,:,i)','EdgeColor','none')

figure,
h = bar(tarDom{f},histc(features.data(inSync(1):inSync(2),f),tarDom{f}),'histc');
h.FaceColor = 'c';
h.EdgeColor = 'c';
h.FaceAlpha = 0.3;
h.EdgeAlpha = 0.3;
hold on
hr = bar(tarDom{f},histc(rfet.data(inSync(1):inSync(2),f),tarDom{f}),'histc');
hr.FaceColor = 'r';
hr.EdgeColor = 'r';
hr.FaceAlpha = 0.3;
hr.EdgeAlpha = 0.3;



eds = linspace(0,140,100);
f= 9;
figure,
h = bar(eds,histc(fet(nniz(fet),f),eds),'histc');
h.FaceColor = 'c';
h.EdgeColor = 'c';
h.FaceAlpha = 0.3;
h.EdgeAlpha = 0.3;
hold on
hr = bar(eds,histc(rfet(nniz(rfet),f),eds),'histc');
hr.FaceColor = 'r';
hr.EdgeColor = 'r';
hr.FaceAlpha = 0.3;
hr.EdgeAlpha = 0.3;
