% req20170914 ----------------------------------------------------
%
%  Status: active
%  Type: Analysis 
%  Final_Forms: NA
%  Description: revisiting repair rigidbody partial reconstruction 
%  Bugs: NA


Session = MTASession.validate('jg05-20120317.cof.all');
Session = MTASession.validate('Ed05-20140529.ont.all');
varargin = {};



markerTriad.nck
nck = 1;
sq(rxyz(ind,reconstructedSolutionsMarkerInds(nck,1),:)-rxyz(ind,markerTriad.nck(nck,2),:))
mtBasis{ind,nck}'*reconstructedSolution.mean{nck,m}'

figure,plot(mrdist(:,:,1,2))

%size(markerTriad.oriMarkerNckCoor) = [522990           4           1           3];
%size(markerTriad.solutions)        = [522990           4           1           3];

nck = 1;
sqrt(sum((markerTriad.solutions(ind,nck,1,:)+sign(s).*markerTriad.oriMarkerNckCoor(ind,nck,1,:)).^2,4))

rhm = fet_rhm(Session);
ncp = fet_ncp(Session);
ncp.unity;

figure,
hold on,
plot(mrdist(:,:,1,2))
plot(rhm.data*2+86);
plot(ncp.data+75);

plot(rhm.data*2+35);
plot(ncp.data+25);

ind = 1:size(rhm,1);
for nck = 1:4,
[th(:,nck),phi(:,nck),r(:,nck)] = ...
    cart2sph(markerTriad.solutions(ind,nck,1,1)+sign(s).*markerTriad.oriMarkerNckCoor(ind,nck,1,1),...
             markerTriad.solutions(ind,nck,1,2)+sign(s).*markerTriad.oriMarkerNckCoor(ind,nck,1,3),...
             markerTriad.solutions(ind,nck,1,3)+sign(s).*markerTriad.oriMarkerNckCoor(ind,nck,1,2));
end


fphi = rhm.copy;
fphi.data = phi;
fphi.filter('ButFilter',3,[4,14],'bandpass');

figure();hold('on');
%plot(nunity(fphi(:,1)).*7)
%plot(rhm.data*4+2);
plot(ncp.data/2-2);
%plot(rhm.data*4+nunity(fphi(:,1)).*7)
ang = create(MTADang,Session,xyz);
plot(ang(:,'head_back','head_front',2)*10);