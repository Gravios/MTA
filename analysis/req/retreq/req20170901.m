% req20170901 ----------------------------------------------------
%
%  Status: inactive
%  Type: Analysis 
%  Final_Forms: NA
%  Description: Find head only features for grooming detection
%  Bugs: NA


Trial = MTATrial('jg05-20120317');

stc = Trial.load('stc','hand_labeled_rev3_jg');


xyz = preproc_xyz(Trial);
% add coordinates of the model's center of mass to the xyz object
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(xyz(:,'hcom',:),3,[5]./(xyz.sampleRate/2),'low'));
xyz.addMarker('flhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(xyz(:,'hcom',:),3,[2]./(xyz.sampleRate/2),'low'));
xyz.addMarker('fvlhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(xyz(:,'hcom',:),3,[0.2]./(xyz.sampleRate/2),'low'));
nz = -cross(xyz(:,'head_back',:)-xyz(:,'hcom',:),xyz(:,'head_left',:)-xyz(:,'hcom',:));
nz = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3)));
nm = nz.*20+xyz(:,'hcom',:);
xyz.addMarker('htx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);

ang = create(MTADang,Trial,xyz);

rhm = fet_rhm(Trial);





hroll = fet_roll(Trial);

hfet = fet_href(Trial);
hfet = fet_head(Trial);

vxy = xyz.vel([],[1,2]);
vxyz = xyz.vel([],[]);

hfig = figure();
subplot2(4,1,[1:3],1);
hold('on');
plot(nunity(diff(xyz(:,'fhcom',3))));
% $$$ plot(nunity(ang(:,'head_back','head_front',2)));
%plot(nunity(sqrt(sum(diff(xyz(:,{'flhcom','head_right'},[1,2]),1,2).^2,3))));
plot(nunity(circ_dist(ang(:,'fhcom','head_right',1),ang(:,'fhcom','head_front',1))));
%plot(nunity(hroll.data));
plot(nunity(vxy(:,'fvlhcom')));
%plot(nunity(hfet(:,[1,2,5])));

subplot2(4,1,[4],1);
plotSTC(stc); 
linkaxes(findobj(hfig,'Type','axes'),'x');



plot(nunity(diff(xyz(:,'hcom',3))));
plot(nunity(sqrt(sum(diff(xyz(:,{'fhcom','head_right'},[1,2]),1,2).^2,3))));







plot(nunity(vxy(:,'fvlhcom')));
tfet = [nunity([0;diff(xyz(:,'hcom',3))]),...
        nunity(circ_dist(ang(:,'fhcom','head_right',1),ang(:,'fhcom','head_front',1))),...
                        nunity(sqrt(sum(diff(xyz(:,{'fhcom','head_right'},[1,2]),1,2).^2,3)))];
wfet = ones([size(xyz,1),size(tfet,2)]);
wfet(nniz(tfet),:) = WhitenSignal(tfet(nniz(tfet),:));

[ys,fs,ts] = mtchglong(wfet,...
                       2^8,...
                       xyz.sampleRate,...
                       2^7,...
                       2^7*0.875,[],[],[],[2,20]);

hfig = figure();
subplot(411);
imagesc(ts,fs,log10(ys(:,:,1,1)'));
caxis([-5,-1]);
ylim([2,20]);
axis xy;
subplot(412);
imagesc(ts,fs,log10(ys(:,:,2,2)'));
caxis([-5.5,-3]);
ylim([2,20]);
axis xy;
subplot(413);
plot([1:size(vxy)]./xyz.sampleRate,vxy(:,'fvlhcom'));
subplot(414);
plotSTC(stc,1);
linkaxes(findobj(hfig,'Type','axes'),'x');