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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Trial = MTATrial.validate('Ed03-20140624'); goodIndex = 16250;goodIndex = 16250;
% LOAD xyz data
xyz = Trial.load('xyz');

% COMPUTE center of mass of the head
headRigidBodyMarkers={'head_back','head_left','head_front','head_right','head_top'};
headRigidBodyMarkers={'head_back','head_left','head_front','head_right'};
headRigidBody = Trial.xyz.model.rb(headRigidBodyMarkers);
hcom = xyz.com(headRigidBody);
hxyz = xyz.copy;
hxyz.model = headRigidBody;
hxyz.data = xyz(:,headRigidBodyMarkers,:);

% ADD center of mass marker to xyz object and a lowpass copy
xyz.addMarker('fhcom',[128,255,128],{{'head_back','head_front',[0,0,1]}},...
               ButFilter(hcom,3,[2]./(Trial.xyz.sampleRate/2),'low'));
xyz.addMarker('hcom',[128,255,128],{{'head_back','head_front',[0,0,1]}},hcom);


% get 
markerIndNCK = nchoosek(1:size(hxyz,2),3);
markerTrioCOM = nan([size(hxyz,1),1,size(hxyz,3)]);
markerTrioCRS = nan([size(hxyz,1),size(markerIndNCK,1),size(hxyz,3)]);
for nck = 1:size(markerIndNCK,1),
    markerTrioCRS(:,nck,:) = cross(hxyz(:,markerIndNCK(nck,1),:)-hxyz(:,markerIndNCK(nck,2),:),...
                                   hxyz(:,markerIndNCK(nck,3),:)-hxyz(:,markerIndNCK(nck,2),:));
    markerTrioCOM(:,nck,:) = mean(hxyz(:,markerIndNCK(nck,:),:),2);

end

markerTrioANG = nan([size(hxyz,1),size(markerIndNCK,1).^2-size(markerIndNCK,1)]);
k = 1;
for i = 1:size(markerIndNCK,1)-1,
    for j = i+1:size(markerIndNCK,1),    
        markerTrioANG(:,k) = acos(dot(markerTrioCRS(:,i,:),markerTrioCRS(:,j,:),3)./...
                                      prod(sqrt(sum(markerTrioCRS(:,[i,j],:).^2,3)),2)).*...
                                      sign(prod(markerTrioCRS(:,[i,j],3),2));
    k = k+1;
    end
end

figure,plot(markerTrioANG(:,:))

% label angular deviation with expectation-maximization of a gaussian hidden markov model 
nind = nniz(markerTrioANG(:,1));
[emgmLabel,emgmModel,emgmLlh] = mixGaussEm(markerTrioANG(nind,[1,2,3])',20);

trgt = markerTrioANG(nind,:);
ul = unique(emgmLabel);
c = 'rgbm';
c = jet(numel(ul));
figure,hold on
for i = ul
    scatter(trgt(emgmLabel'==i,1),trgt(emgmLabel'==i,2),5,c(i,:))
end

emLabel = nan([size(markerTrioANG,1),1]);
emLabel(nind)=emgmLabel



figure,plot(sqrt(sum((markerTrioCOM-circshift(markerTrioCOM,1)).^2,3)))

comDistMat = diff(sqrt(sum((markerTrioCOM-circshift(markerTrioCOM,1)).^2,3)));



% $$$ 
% $$$  
% $$$ nz = -cross(xyz(:,'head_back',:)-hcom,xyz(:,'head_left',:)-hcom);
% $$$ nz = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3)));
% $$$ nm = nz.*20+hcom;
% $$$ xyz.addMarker('htx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);
% $$$ 
% $$$ ny = cross(xyz(:,'htx',:)-hcom,xyz(:,'head_back',:)-hcom);
% $$$ ny = bsxfun(@rdivide,ny,sqrt(sum((ny).^2,3)));
% $$$ nm = ny.*20+hcom;
% $$$ xyz.addMarker('hrx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);
% $$$ 
% $$$ nx = cross(xyz(:,'hrx',:)-hcom,xyz(:,'htx',:)-hcom);
% $$$ nx = bsxfun(@rdivide,nx,sqrt(sum((nx).^2,3)));
% $$$ nm = nx.*20+hcom;    
% $$$ xyz.addMarker('hbx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);
% $$$ 
% $$$ % create surrogate markers before this step !!!
% $$$ hrb = xyz(:,{'chff','chrf','chtf'},:)-xyz(:,{'chcf'},:);
% $$$ hrb = mat2cell(hrb,ones([size(xyz,1),1]),ones([3,3]));
% $$$ U = cell(size(hrb));
% $$$ S = cell(size(hrb));
% $$$ V = cell(size(hrb));
% $$$ [U,S,V] = cellfun(@svd,hrb);
% $$$ V = cell2mat(V);
% $$$ 
% $$$ hyz = sq(xyz(:,'head_back',:))*cell2mat(V)
% $$$ hyz = multiprod(sq(xyz(:,'chf',:)-xyz(:,'chcf',:)),V(:,[2,3],:)',[1,2],2);
% $$$ hfa = atan2(hyz(:,1),hyz(:,2));
% $$$ 
% $$$ 
% $$$ figure,plot(hfa)
% $$$ 
% $$$ 
% $$$ c = sqrt(sum(cross(xyz(:,6,:)-xyz(:,5,:),xyz(:,7,:)-xyz(:,5,:)).^2,3));
% $$$ figure,plot((c-mean(c(nniz(c(:)))))./std(c(nniz(c(:)))))
% $$$ 
% $$$ 
% $$$ b = sqrt(sum(cross(xyz(:,6,:)-xyz(:,8,:),xyz(:,7,:)-xyz(:,8,:)).^2,3));
% $$$ 
% $$$ 

% $$$ 




Trial = MTATrial.validate('jg05-20120317.cof.all');
xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,4,'low');
ang = create(MTADang,Trial,xyz);

dz = MTADxyz('data',log10(1+abs([0;diff(xyz(:,'head_front',3))])),'sampleRate',xyz.sampleRate);

figure
edx = linspace(0,.8,50);
edy = linspace(-1,pi/2,50);
ind = Trial.stc{'a'};
hist2([dz(ind),ang(ind,'spine_middle','spine_upper',2)],edx,edy);
caxis([0,500])
figure
ind = Trial.stc{'a-r'};
hist2([dz(ind),ang(ind,'spine_middle','spine_upper',2)],edx,edy);
caxis([0,500])

figure,
ind = Trial.stc{'r'};
hist2([dz(ind),ang(ind,'spine_middle','spine_upper',2)],50,50);
caxis([0,500])