Trial = MTATrial('jg05-20120317');

xyz = Trial.load('xyz');

% create a ridgid body model
rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
% find the center of mass of the model
hcom = xyz.com(rb);
% add coordinates of the model's center of mass to the xyz object
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(hcom,3,[2.4]./(Trial.xyz.sampleRate/2),'low'));
xyz.addMarker('lfhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              circshift(ButFilter(hcom,3,[.8]./(Trial.xyz.sampleRate/2),'low'),30));
xyz.addMarker('vlfhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              circshift(ButFilter(hcom,3,[.2]./(Trial.xyz.sampleRate/2),'low'),120));
xyz.addMarker('vlfhcomn',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(hcom,3,[.2]./(Trial.xyz.sampleRate/2),'low'));

xyz.filter('ButFilter',3,30,'low');


ang = create(MTADang,Trial,xyz);

vxy = xyz.vel([],[1,2]);



cang = circ_dist(ang(:,'fhcom','lfhcom',1),ang(:,'head_front','lfhcom',1));
dcang = circ_dist(circshift(cang,10),circshift(cang,-10));


bang = MTADxyz('data',clip(ang(:,'fhcom','lfhcom',3)./sqrt(sq(sum(GetSegs(circshift(dcang,15),1:size(dcang),30).^2)))'./10,0,1000),'sampleRate',ang.sampleRate);
bang.filter('ButFilter',3,[2.4],'low');

[thl,phl,rhl] = cart2sph(xyz(:,'vlfhcom',1)-circshift(xyz(:,'vlfhcom',1),0),xyz(:,'lfhcom',2)-circshift(xyz(:,'vlfhcom',2),0),xyz(:,'lfhcom',3)-circshift(xyz(:,'vlfhcom',3),0));

[thl,phl,rhf] = cart2sph(xyz(:,'lfhcom',1)-circshift(xyz(:,'lfhcom',1),1),xyz(:,'lfhcom',2)-circshift(xyz(:,'lfhcom',2),1),xyz(:,'lfhcom',3)-circshift(xyz(:,'lfhcom',3),1));


m = 'lfhcom';
dtan = [xyz(:,m,1)-circshift(xyz(:,m,1),-round(0.1*xyz.sampleRate)),...
        xyz(:,m,2)-circshift(xyz(:,m,2),-round(0.1*xyz.sampleRate)),...
        xyz(:,m,3)-circshift(xyz(:,m,3),-round(0.1*xyz.sampleRate))];
%[tanTh,tanPh,tanRh] = cart2sph(dtan(:,1),dtan(:,2),dtan(:,3));

m = 'head_front';
dhed = [xyz(:,m,1)-xyz(:,'head_back',1),...
        xyz(:,m,2)-xyz(:,'head_back',2),...
        xyz(:,m,3)-xyz(:,'head_back',3)];


bpr = MTADxyz('data',log10(clip(-sum(dtan.*bsxfun(@rdivide,dhed,sqrt(sum(dhed.^2,2))),2),-19,200)+20),'sampleRate',xyz.sampleRate);


eds = linspace([1.2,2.3]);
figure,hold on
ind = Trial.stc{'a-w'};
ha = bar(eds,histc(bpr(ind),eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .4;
ha.EdgeColor = 'c';
ha.EdgeAlpha = .4;
ind = Trial.stc{'w'};
hs = bar(eds,histc(bpr(ind),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeColor = 'r';
hs.EdgeAlpha = .4;



ind = 1000:2000;
figure,
plot(xyz(ind,'vlfhcom',1),xyz(ind,'vlfhcom',2))
hold on
plot(xyz(ind,'lfhcom',1),xyz(ind,'lfhcom',2))
hold on
plot(xyz(ind,'fhcom',1),xyz(ind,'fhcom',2))


figure,plot(bang.data)
hold on,plot(vxy(:,1))
hold on,plot(vxy(:,'head_front'))
hold on,plot(rh.*100);
hold on,plot(ang(:,'vlfhcom','lfhcom',3))
Lines(Trial.stc{'n'}(:),[],'g');
Lines(Trial.stc{'w'}(:),[],'b');


lbang = bang.copy;
lbang.data(lbang(:,1)<1e-6) = 1e-6;
lbang.data = log10(lbang.data);

eds = linspace([-6,3,200]);
figure,hold on
ind = Trial.stc{'a-r-m'};
ha = bar(eds,histc(lbang(ind),eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .4;
ha.EdgeColor = 'c';
ha.EdgeAlpha = .4;
ind = Trial.stc{'m'};
hs = bar(eds,histc(lbang(ind),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeColor = 'r';
hs.EdgeAlpha = .4;



fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.4,'low');
lvxy = fxyz.vel([],[1,2]);
lvxy.data(vxy(:,1)<1e-3) = 1e-3;
lvxy.data = log10(lvxy.data);

eds = linspace([-3,3,100]);
figure,hold on
ind = Trial.stc{'a-w-r'};
ha = bar(eds,histc(lvxy(ind,1),eds),'histc');
ha.FaceColor = 'c';
ha.FaceAlpha = .4;
ha.EdgeColor = 'c';
ha.EdgeAlpha = .4;
ind = Trial.stc{'w'};
hs = bar(eds,histc(lvxy(ind,1),eds),'histc');
hs.FaceColor = 'r';
hs.FaceAlpha = .4;
hs.EdgeColor = 'r';
hs.EdgeAlpha = .4;


