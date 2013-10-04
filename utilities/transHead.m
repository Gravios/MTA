
sl
hr = s.xyz(:,s.Model.gmi('head_right'),:)-s.xyz(:,s.Model.gmi('head_back'),:);
hl = s.xyz(:,s.Model.gmi('head_left'),:)-s.xyz(:,s.Model.gmi('head_back'),:);
crl = sq(cross(hl,hr));
crlt = reshape(repmat(sqrt(sum(crl.^2,2)),3,1),[],3);
trans_xyz = crl./crlt;
hcom = s.com(s.Model.rb({'head_back','head_left','head_front','head_right'}));

%xyz = hcom+reshape(trans_xyz,[],1,3).*20;
%xyz = s.xyz(:,s.Model.gmi('head_front'),:)+reshape(trans_xyz,[],1,3).*20;
%xyz = s.xyz(:,s.Model.gmi('head_back'),:)+reshape(trans_xyz,[],1,3).*20;
s.xyz(:,s.Model.gmi({'head_back','head_left','head_front','head_right'}),:) = s.xyz(:,s.Model.gmi({'head_back','head_left','head_front','head_right'}),:)+reshape(repmat(trans_xyz,4,1),[],4,3).*20;
s.xyz(isnan(s.xyz)) = 0;

name = 'head_front_td';
color = '[0,0,1]';
sticks = {{'spine_upper',name,[0,0,1]},{'head_back',name,[0,0,1]},{'head_left',name,[0,1,0]},{'head_front',name,[0.5,0,0.5]},{'head_right',name,[1,0,0]}};
s = s.addMarker(name,color,sticks,xyz);

