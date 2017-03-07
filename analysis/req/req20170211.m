% NAME :    req20170211
% PURPOSE : play with respiration timeseries reconstruction with nn
%           
%


Trial = MTATrial.validate('Ed05-20140529.ont.all');
ncp = fet_ncp(Trial,[],[],2);
ncp.filter('ButFilter',3,50,'low');

ncp.data = clip(ncp.data,-8000,8000);

eds = linspace(-8000,8000,200);
ncpb = discretize(ncp.data,eds);

rhm = fet_rhm(Trial);

xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,30,'low');
ang = create(MTADang,Trial,xyz);

fet = [];
fet = [fet,sq(ang(:,1,2,:))];
fet = [fet,sq(ang(:,1,3,:))];
fet = [fet,sq(ang(:,1,4,:))];
fet = [fet,sq(ang(:,1,7,:))];
fet = [fet,sq(ang(:,2,3,:))];
fet = [fet,sq(ang(:,2,4,:))];
fet = [fet,sq(ang(:,2,7,:))];
fet = [fet,sq(ang(:,3,4,:))];
fet = [fet,sq(ang(:,3,7,:))];
fet = [fet,sq(ang(:,4,7,:))];
fet = [fet,rhm.data];

nind = nniz(fet);

net = patternnet(200);
%net = feedforwardnet(1000);

ncpmat = zeros([numel(eds),size(ncp,1)]);
ncpmat((0:size(ncp,1)-1)'*numel(eds)+ncpb)=1;
ncpmat = ncpmat';

[net,tr] = train(net,fet(nind,:)',ncpmat(nind,:)');

out = zeros([size(ncp,1),numel(eds)]);
out(nind,:) = net(fet(nind,:)')';