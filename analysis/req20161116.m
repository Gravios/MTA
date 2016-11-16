Trial = MTATrial.validate('jg05-20120317');
Trial = MTATrial.validate('Ed03-20140624');

xyz = Trial.load('xyz');
hcom = xyz.com(Trial.xyz.model.rb({'head_back','head_left','head_front','head_right'}));
xyz.addMarker('fhcom',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(hcom,3,[2]./(Trial.xyz.sampleRate/2),'low'));
xyz.addMarker('fhf',[.7,1,.7],{{'head_back','head_front',[0,0,1]}},ButFilter(xyz(:,7,:),3,[2]./(Trial.xyz.sampleRate/2),'low'));

ang = create(MTADang,Trial,xyz);


ncp = fet_ncp(Trial);

figure,
plot(ang(:,'fhcom','head_front',2)-ang(:,'fhcom','fhf',2)),
hold on,
plot(ang(:,'fhcom','head_front',1)-ang(:,'fhcom','fhf',1)),
plot((ang(:,'fhcom','head_front',3)-mean(ang(nniz(ang(:,'fhcom','head_front',3)),'fhcom','head_front',3)))./10)
plot(nunity(ncp.data)./5)

figure,hold on,
plot(diff(ButFilter(ang(:,'fhcom','head_front',2)-ang(:,'fhcom','fhf',2),3,[1,20]./(Trial.ang.sampleRate/2),'bandpass')))
plot(diff(ButFilter(ang(:,'fhcom','head_front',1)-ang(:,'fhcom','fhf',1),3,[1,20]./(Trial.ang.sampleRate/2),'bandpass')))
plot(diff(ButFilter((ang(:,'fhcom','head_front',3)-mean(ang(nniz(ang(:,'fhcom','head_front',3)),'fhcom','head_front',3)))./10,3,[1,20]./(Trial.ang.sampleRate/2),'bandpass')))
plot(nunity(ncp.data)./50)



spitch = diff(ButFilter(ang(:,'fhcom','head_front',2)-ang(:,'fhcom','fhf',2),3,[1,20]./(Trial.ang.sampleRate/2),'bandpass'));
sasmth = diff(ButFilter(ang(:,'fhcom','head_front',1)-ang(:,'fhcom','fhf',1),3,[1,20]./(Trial.ang.sampleRate/2),'bandpass'));
rhm = diff(ButFilter(ang(:,'fhcom','head_front',3),3,[1,20]./(Trial.ang.sampleRate/2),'bandpass'));
mncp = nunity(ncp.data(1:end-1))./50;


[ys,fs,ts] = fet_spec(Trial,MTADxyz('data',[spitch,sasmth,rhm,mncp],'sampleRate',xyz.sampleRate),'mtchglong',false);

figure,
for i = 1:size(ys,3),
    sp(i)=subplot(size(ys,3),1,i);
    imagesc(ts,fs,ys(:,:,i,4)');    
    axis xy
end
linkaxes(sp,'xy');
caxis([0,1e-6])



figure,
sp(1)=subplot(211); imagesc(ts,fs,log10(ys(:,:,3,3))');axis xy;caxis([-4,-2.2])
sp(2)=subplot(212); imagesc(ts,fs,log10(ys(:,:,4,4))');axis xy;caxis([-9,-3.8])
linkaxes(sp,'xy');



% create surrogate markers before this step !!!
hrb = xyz(:,{'chff','chrf','chtf'},:)-xyz(:,{'chcf'},:);
hrb = mat2cell(hrb,ones([size(xyz,1),1]),ones([3,3]));
U = cell(size(hrb));
S = cell(size(hrb));
V = cell(size(hrb));
[U,S,V] = cellfun(@svd,hrb);
V = cell2mat(V);

hyz = sq(xyz(:,'head_back',:))*cell2mat(V)
hyz = multiprod(sq(xyz(:,'chf',:)-xyz(:,'chcf',:)),V(:,[2,3],:)',[1,2],2);
hfa = atan2(hyz(:,1),hyz(:,2));


figure,plot(hfa)


c = sqrt(sum(cross(xyz(:,6,:)-xyz(:,5,:),xyz(:,7,:)-xyz(:,5,:)).^2,3));
figure,plot((c-mean(c(nniz(c(:)))))./std(c(nniz(c(:)))))


b = sqrt(sum(cross(xyz(:,6,:)-xyz(:,8,:),xyz(:,7,:)-xyz(:,8,:)).^2,3));






