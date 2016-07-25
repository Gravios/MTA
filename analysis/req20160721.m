Trial = MTATrial('jg05-20120317');

xyz = Trial.load('xyz');

% create a ridgid body model
rb = Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'});
% find the center of mass of the model
hcom = xyz.com(rb);
% add coordinates of the model's center of mass to the xyz object
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(hcom,3,[2.4]./(Trial.xyz.sampleRate/2),'low'));
xyz.addMarker('lfhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              circshift(ButFilter(hcom,3,[.5]./(Trial.xyz.sampleRate/2),'low'),30));
xyz.filter('ButFilter',3,30,'low');


ang = create(MTADang,Trial,xyz);

fet = Trial.xyz.copy;

vxy = xyz.vel([],[1,2]);

bang = [ButFilter(ang(:,'head_back','fhcom',3) ,3,[1,50]./(Trial.ang.sampleRate/2),'bandpass'),...
        ButFilter(ang(:,'head_right','fhcom',3),3,[1,50]./(Trial.ang.sampleRate/2),'bandpass')];
%ButFilter(ang(:,'head_top','fhcom',3),3,[2,50]./(Trial.ang.sampleRate/2),'bandpass')];


pang = [circ_dist(circshift(ang(:,'head_back','fhcom',1),0),circshift(ang(:,'head_back','fhcom',1),-1)),...
        circ_dist(circshift(ang(:,'head_right','fhcom',1),0),circshift(ang(:,'head_right','fhcom',1),-1))];

fet.data = [0,0;ButFilter(diff(bang),3,[2,50]/(ang.sampleRate/2),'bandpass')];

figure,plot(fet.data)
figure,plot(bang)
%figure,plot(pang)
Lines(Trial.stc{'n'}(:),[],'g');
Lines(Trial.stc{'w'}(:),[],'b');


cang = circ_dist(ang(:,'fhcom','lfhcom',1),ang(:,'head_front','lfhcom',1));
dcang = circ_dist(circshift(cang,10),circshift(cang,-10));

figure, plot(ang(:,'fhcom','lfhcom',3))
hold on,plot(vxy(:,1))
hold on,plot(dcang.*40)

figure,plot(clip(ang(:,'fhcom','lfhcom',3)./sqrt(sq(sum(GetSegs(circshift(dcang,15),1:size(dcang),30).^2)))'./10,0,1000))
hold on,plot(vxy(:,1))
hold on,plot(ang(:,'fhcom','lfhcom',3))
Lines(Trial.stc{'n'}(:),[],'g');
Lines(Trial.stc{'w'}(:),[],'b');
