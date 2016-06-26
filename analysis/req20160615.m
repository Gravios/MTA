Trial = MTATrial.validate('Ed03-20140624.cof.all');


xyzh = [ 23.71, 14.43, 74.78;... head_back
         52.44,-11.00, 66.56;... head_left
         76.82,  4.47, 65.67;... head_front
         51.46, 35.80, 65.67;... head_right
         79.47, 10.00, 31.12 ... head_nose
];

cxyzh = mean(xyzh(1:3,:));

xyzh = bsxfun(@minus,xyzh,cxyzh);

xyz = Trial.load('xyz');

xyz.addMarker('hcom',...     Name
              [.7,0,.7],...  Color
              {{'head_back', 'hcom',[0,0,255]},... Sticks to visually connect
               {'head_left', 'hcom',[0,0,255]},... new marker to skeleton
               {'head_front','hcom',[0,0,255]},...
               {'head_right','hcom',[0,0,255]}},... 
              xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));

ang = create(MTADang,Trial,xyz);
% $$$ 
% $$$ markers = {'head_back','head_left','head_front','head_right'};              
% $$$ A = xyz(:,markers,:);
% $$$ A = bsxfun(@minus,A,mean(A));
% $$$ A = permute(A,[2,3,1]);
% $$$ A = sq(mat2cell(A,size(A,1),size(A,2),ones([size(A,3),1])));
% $$$ [U,S,V]=cellfun(@svd,A,'uniformoutput',false);
% $$$ 
% $$$ 
% $$$ nmar =cellfun(@mtimes,V,repmat({[30,0,20]'},size(A)),'uniformoutput',false);
% $$$ 
% $$$ xyz.addMarker('hb_syn',...     Name
% $$$               [.7,0,.7],...  Color
% $$$               {{'head_back', 'hcom',[0,0,255]},... Sticks to visually connect
% $$$                {'head_left', 'hcom',[0,0,255]},... new marker to skeleton
% $$$                {'head_front','hcom',[0,0,255]},...
% $$$                {'head_right','hcom',[0,0,255]}},... 
% $$$               permute(reshape(cell2mat(nmar),3,[]),[2,3,1])+xyz(:,'hcom',:));

              
nz = cross(xyz(:,'head_front',:)-xyz(:,'hcom',:),xyz(:,'head_left',:)-xyz(:,'hcom',:));
nz = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3)));
nm = nz.*20+xyz(:,'hcom',:);
xyz.addMarker('htx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);

ny = cross(xyz(:,'htx',:)-xyz(:,'hcom',:),xyz(:,'head_front',:)-xyz(:,'hcom',:));
ny = bsxfun(@rdivide,ny,sqrt(sum((ny).^2,3)));
nm = ny.*20+xyz(:,'hcom',:);
xyz.addMarker('hrx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);

nx = cross(xyz(:,'hrx',:)-xyz(:,'hcom',:),xyz(:,'htx',:)-xyz(:,'hcom',:));
nx = bsxfun(@rdivide,nx,sqrt(sum((nx).^2,3)));
nm = nx.*20+xyz(:,'hcom',:);    
xyz.addMarker('hbx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);

nc = cat(2,nx,ny,nz);
A = permute(nc,[3,2,1]);
A = sq(mat2cell(A,size(A,1),size(A,2),ones([size(A,3),1])));

nmar =cellfun(@mtimes,A,repmat({[-30,0,-20]'},size(A)),'uniformoutput',false);

xyz.addMarker('hb_nc',...     Name
              [.7,0,.7],...  Color
              {{'head_back', 'hcom',[0,0,255]},... Sticks to visually connect
               {'head_left', 'hcom',[0,0,255]},... new marker to skeleton
               {'head_front','hcom',[0,0,255]},...
               {'head_right','hcom',[0,0,255]}},... 
              permute(reshape(cell2mat(nmar),3,[]),[2,3,1])+xyz(:,'hcom',:));

              
              
              
oxyz = Trial.load('xyz');
oxyz.addMarker('hcom',...     Name
              [.7,0,.7],...  Color
              {{'head_back', 'hcom',[0,0,255]},... Sticks to visually connect
               {'head_left', 'hcom',[0,0,255]},... new marker to skeleton
               {'head_front','hcom',[0,0,255]},...
               {'head_right','hcom',[0,0,255]}},... 
              oxyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));
foxyz = oxyz.copy;
foxyz.filter('ButFilter',3,2,'low');
oxyz.addMarker('fhcom',...     Name
              [.7,0,.7],...  Color
              {{'head_back', 'fhcom',[0,0,255]},... Sticks to visually connect
               {'head_left', 'fhcom',[0,0,255]},... new marker to skeleton
               {'head_front','fhcom',[0,0,255]},...
               {'head_right','fhcom',[0,0,255]}},... 
              foxyz(:,'hcom',:));
              
txyz = Trial.load('xyz','trb');
txyz.addMarker('hcom',...     Name
              [.7,0,.7],...  Color
              {{'head_back', 'hcom',[0,0,255]},... Sticks to visually connect
               {'head_left', 'hcom',[0,0,255]},... new marker to skeleton
               {'head_front','hcom',[0,0,255]},...
               {'head_right','hcom',[0,0,255]}},... 
              txyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));
ftxyz = txyz.copy;
ftxyz.filter('ButFilter',3,2,'low');
txyz.addMarker('fhcom',...     Name
              [.7,0,.7],...  Color
              {{'head_back', 'fhcom',[0,0,255]},... Sticks to visually connect
               {'head_left', 'fhcom',[0,0,255]},... new marker to skeleton
               {'head_front','fhcom',[0,0,255]},...
               {'head_right','fhcom',[0,0,255]}},... 
              ftxyz(:,'hcom',:));

oang = create(MTADang,Trial,oxyz);
tang = create(MTADang,Trial,txyz);

ncp = fet_ncp(Trial,oxyz);
rhm = fet_rhm(Trial,oxyz);
trhm = fet_rhm(Trial,oxyz);

figure,
plot(oang(:,5,'fhcom',3))
hold on
plot(tang(:,5,'fhcom',3))
plot(nunity(ncp.data)+24)
plot(nunity(rhm.data)*3+24)
plot(nunity(trhm.data)*3+24)
plot(tang(:,5,7,2)*10+28)