Trial = MTATrial('jg05-20120317');
Trial = MTATrial('jg05-20120317');
Trial = MTATrial('jg03-20110501');

xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,2.5,'low');
ang = create(MTADang,Trial,xyz);

xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,.5,'low');
ang = create(MTADang,Trial,xyz);

figure,
sp = []
sp(1) = subplot2(4,1,1:3,1);
plot(circ_dist(circshift(ang(:,1,4,1),1),circshift(ang(:,1,4,1),-1)));
grid on
sp(2) = subplot2(4,1,4,1);
plotSTC(Trial.stc,ang.sampleRate,[],{'walk','rear','turn','pause'},'wrgc');
grid on
linkaxes(sp,'x');

Trial.load('stc','hand_labeled_rev3_jg')


cpang = [];
cdadt = [];
figure, hold on
for turn = Trial.stc{'n',ang.sampleRate}.data'
    pang = ang(turn',1,4,1);
    pang = circ_dist(pang,pang(1));
    dadt = circ_dist(circshift(pang,1),circshift(pang,-1));
    plot(abs(pang),abs(dadt),'.');
    cpang = [cpang;pang];
    cdadt = [cdadt;dadt];
end

figure
hist2([cpang,log10(cdadt)],linspace([0,pi,30]),linspace([-3,-1,30]))

figure, hold on
for turn = Trial.stc{'w',ang.sampleRate}.data'
    pang = ang(turn',1,4,1);
    pang = circ_dist(pang,pang(1));
    dadt = circ_dist(circshift(pang,1),circshift(pang,-1));
    plot(abs(pang),abs(dadt),'.');
    cpang = [cpang;pang];
    cdadt = [cdadt;dadt];
end

figure
hist2([cpang,log10(cdadt)],linspace([0,pi,30]),linspace([-3,-1,30]))

figure,
subplot(211),hist(diff(Trial.stc{'w',1}.data,1,2),30)
subplot(212),hist(diff(Trial.stc{'n',1}.data,1,2),30)


Trial = MTATrial('jg05-20120317');
xyz = Trial.load('xyz');
rb = xyz.model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'});
% find the center of mass of the model
hcom = xyz.com(rb);
% add coordinates of the model's center of mass to the xyz object
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(hcom,4,[1]./(xyz.sampleRate/2),'low'));
ang = create(MTADang,Trial,xyz);
Trial.load('stc','hand_labeled_rev3_jg')


figure,hold on
plot(ang(:,4,10,2))
plot(ang(:,3,10,2))
plot(ang(:,2,10,2))
plot(ang(:,1,10,2))
Lines(Trial.stc{'w'}(:),[],'k');
Lines(Trial.stc{'r'}(:),[],'r');

figure,hold on
plot(circ_dist(ang(:,1,10,1),ang(:,10,4,1)))
plot(circ_dist(ang(:,2,10,1),ang(:,10,4,1)))
plot(circ_dist(ang(:,2,10,1),ang(:,10,3,1)))
plot(circ_dist(ang(:,1,10,1),ang(:,10,3,1)))
Lines(Trial.stc{'w'}(:)-.5,[],'k');
Lines(Trial.stc{'r'}(:)+.5,[],'r');
Lines(Trial.stc{'n'}(:),[],'g');




swag = MTADang('data',[circ_dist(ang(:,1,10,1),ang(:,10,4,1)),...
                       circ_dist(ang(:,2,10,1),ang(:,10,4,1)),...
                       circ_dist(ang(:,2,10,1),ang(:,10,3,1)),...
                       circ_dist(ang(:,1,10,1),ang(:,10,3,1))],...
               'sampleRate',ang.sampleRate);

swag.data(~nniz(xyz),:) = 0;

swag.filter('ButFilter',3,1,'high');
figure,plot(swag.data)
Lines(Trial.stc{'w'}(:)-.5,[],'k');
Lines(Trial.stc{'r'}(:)+.5,[],'r');
Lines(Trial.stc{'n'}(:),[],'g');


figure,plot(prod(swag.data,2))




%% Walk golani 1980
QuickSessionSetup(get_session_list('jg03'));

s = MTASession('jg03-20110501');

Trial = MTATrial('jg03-20110501');
Trial.load('stc','NN0317');

xyz = Trial.load('xyz');
rb = xyz.model.rb({'spine_lower','pelvis_root','spine_middle','spine_upper'});
hcom = xyz.com(rb);
xyz.addMarker('fbcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(hcom,4,[1]./(xyz.sampleRate/2),'low'));
rb = xyz.model.rb({'head_back','head_left','head_front','head_right'});
hcom = xyz.com(rb);
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},...
              ButFilter(hcom,4,[1]./(xyz.sampleRate/2),'low'));
clear('hcom');

ang = create(MTADang,Trial,xyz);
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,5,'low');
fvxy = fxyz.vel([],[1,2]);
hfvxy = fvxy.copy;
hfvxy.filter('ButFilter',3,1,'high');
fang = create(MTADang,Trial,fxyz);
hxyz = xyz.copy;
hxyz.filter('ButFilter',3,1,'high');
hfxyz = fxyz.copy;
hfxyz.filter('ButFilter',3,1,'high');


figure,
plot(circ_dist(fang(:,'spine_lower','spine_middle',1),...
               fang(:,'hip_right','hip_left',1))...
     );


figure,hold on,
plot(nunity(ang(:,'hip_right','hip_left',2)));
plot(nunity(fang(:,'hip_right','hip_left',2)));
plot(nunity(fxyz(:,'hip_left',3)));
plot(nunity(fxyz(:,'hip_right',3)));
plot(nunity(fxyz(:,'spine_lower',3)));
plot(nunity(fvxy(:,'pelvis_root')));
plot(nunity(fvxy(:,'spine_lower')));
plot(nunity(fvxy(:,'hip_left')),'-.');
plot(nunity(fvxy(:,'hip_right')),'.');
Lines(Trial.stc{'w'}(:)-.5,[],'k');
Lines(Trial.stc{'r'}(:)+.5,[],'r');

figure,hold on,
plot(nunity(hfvxy(:,'spine_lower')),'b')
plot(nunity(hfvxy(:,'hip_left')),'m-.')
plot(nunity(hfvxy(:,'hip_right')),'g.')
plot(nunity(hfvxy(:,'hip_left'))-nunity(hfvxy(:,'hip_right')),'c-.')
plot(nunity(hfvxy(:,'spine_middle')))
plot(nunity(circ_dist(circshift(ang(:,'spine_lower','spine_middle',1),10),circshift(ang(:,'spine_lower','spine_middle',1),-10))),'k')
plot(circ_dist(fang(:,'pelvis_root','fbcom',1),fang(:,'fbcom','spine_middle',1)),'r.')
plot(nunity(hfxyz(:,'spine_lower',3)),'m');
Lines(Trial.stc{'w'}(:)-.5,[],'k');
Lines(Trial.stc{'n'}(:),[],'g');
     
figure, hold on
plot(circ_dist(fang(:,'spine_lower','spine_middle',1),fang(:,'spine_lower','pelvis_root',1)),'-')
plot(circ_dist(fang(:,'spine_lower','spine_middle',1),fang(:,'spine_lower','hip_left',1)),'-.')
plot(circ_dist(-fang(:,'spine_lower','spine_middle',1),fang(:,'spine_lower','hip_right',1)),'-.')
Lines(Trial.stc{'w'}(:)-.5,[],'k');

xlm = [45300,45700];
figure,hold on
subplot(311);hold on
%plot(nunity(hfxyz(:,[2,4],3)));
% $$$ plot(ang(:,1,3,3));
% $$$ plot(ang(:,1,5,3));
plot(ang(:,6,11,2));
plot(fvxy(:,1)); hold on
plot(diff(fxyz(:,1,3)));
plot(ang(:,2,4,2));
xlim(xlm);

% $$$ Lines(Trial.stc{'w'}(:)-.5,[],'k');
% $$$ Lines(Trial.stc{'n'}(:),[],'g');

%plot(nunity(diff(fxyz(:,[2,4],3))));
% $$$ plot(nunity(hfvxy(:,'hip_left')),'m-.')
% $$$ plot(nunity(hfvxy(:,'hip_right')),'g-.')
% $$$ plot(nunity(hfvxy(:,'spine_lower')),'b')
% $$$ plot(nunity(hfvxy(:,'spine_middle')),'r')


tvec = circshift(fxyz(:,1,[1,2]),-10)-circshift(fxyz(:,1,[1,2]),10);
bvec = fxyz(:,4,[1,2])-fxyz(:,1,[1,2]);
ubvec = bsxfun(@rdivide,bvec,sqrt(sum(bvec.^2,3)));
figure,plot(dot(tvec,ubvec,3))

xlim(xlm);

subplot(312);hold on,
plot(circ_dist(fang(:,'pelvis_root','fbcom',1),fang(:,'fbcom','spine_middle',1)),'r.')
plot(circ_dist(fang(:,'spine_lower','fbcom',1),fang(:,'fbcom','spine_upper',1)),'m.')
xlim(xlm);

subplot(313);hold on
plot(nunity(hfxyz(:,[1,3,5,6,7],3)));
legend({'SL','PR','SM','SU','HB'})
xlim(xlm);

figure,hold on
plot(nunity(fang(:,1,10,3)).*10)
plot(nunity(hfxyz(:,1,3)));
xlim([221700,222100]);

figure,hold on,
plot(fxyz(:,2

figure,hold on
plot(ang(:,1,10,3))

m = 6
n = 11;
edx = linspace(-pi/2,pi/2,100)
edy = linspace(50 , 100,100);
figure
subplot(321)
ind =  Trial.stc{'a'};
hist2(sq(ang(ind,n,m,[2,3])),edx,edy)
caxis([0,100])
subplot(322)
ind =  Trial.stc{'r'};
hist2(sq(ang(ind,n,m,[2,3])),edx,edy)
caxis([0,100])
subplot(323)
ind =  Trial.stc{'m'};
hist2(sq(ang(ind,n,m,[2,3])),edx,edy)
caxis([0,100])
subplot(324)
ind =  Trial.stc{'p'};
hist2(sq(ang(ind,n,m,[2,3])),edx,edy)
caxis([0,100])
subplot(325)
ind =  Trial.stc{'n'};
hist2(sq(ang(ind,n,m,[2,3])),edx,edy)
caxis([0,100])
subplot(326)
ind =  Trial.stc{'w'};
hist2(sq(ang(ind,n,m,[2,3])),edx,edy)
caxis([0,100])