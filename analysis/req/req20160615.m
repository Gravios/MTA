Trial = MTATrial.validate('Ed03-20140624.cof.all');


% stereotaxic coordinates of head marker positions
xyzh = [ 23.71, 14.43, 74.78;... head_back
         52.44,-11.00, 66.56;... head_left
         76.82,  4.47, 65.67;... head_front
         51.46, 35.80, 65.67;... head_right
];
% Flipped version
xyzh = [ 23.71, 14.43, 74.78;... head_back
         51.46, 35.80, 65.67;... head_right
         76.82,  4.47, 65.67;... head_front
         52.44,-11.00, 66.56;... head_left         
];

nosexyz = [79.47, 10.00, 31.12];
neckxyz = [79.47-30, 10.00, 31.12];

% Ed05
Trial = MTATrial.validate('Ed05-20140529.ont.all');
xyzh = [ 22.46, 11.96, 61.07;... head_back
         52.42,-18.00, 60.24;... head_left                  
         80.00, -2.38, 57.25;... head_front
         51.48, 34.50, 65.02;... head_right
];


xyzh = [ 22.46, 11.96, 61.07;... head_back
         52.42,-18.00, 60.24;... head_left                  
         80.00, -2.38, 57.25;... head_front
         51.48, 34.50, 65.02;... head_right
];

yshift  = min(xyzh(:,2));
xyzh = bsxfun(@times,...
              bsxfun(@plus,xyzh,[0,abs(yshift),0]),...
              [1,-1,1]);

nosexyz = [85.68, 9.65, 32.39];
neckxyz = [85.68-30, 9.65, 32.39];


nosexyz = ([85.68   , 9.65, 32.39]+[0,abs(yshift),0]).*[1,-1,1];
neckxyz = ([85.68-40, 9.65, 32.39]+[0,abs(yshift),0]).*[1,-1,1];

% corrected
nosexyz = [85.68, 9.65, 82.11];
neckxyz = [85.68-30, 9.65, 82.11];

2*57.25 - 32.39


% LOAD xyz data
xyz       = preproc_xyz(Trial);

nz = -cross(xyz(:,'head_back',:)-hcom,xyz(:,'head_left',:)-hcom);
nz = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3)));
nm = nz.*20+hcom;
xyz.addMarker('htx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);

ny = cross(xyz(:,'htx',:)-hcom,xyz(:,'head_back',:)-hcom);
ny = bsxfun(@rdivide,ny,sqrt(sum((ny).^2,3)));
nm = ny.*20+hcom;
xyz.addMarker('hrx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);

nx = cross(xyz(:,'hrx',:)-hcom,xyz(:,'htx',:)-hcom);
nx = bsxfun(@rdivide,nx,sqrt(sum((nx).^2,3)));
nm = nx.*20+hcom;    
xyz.addMarker('hbx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);


rigidBody = xyz.model.rb({'head_back','head_left','head_front','head_right'});

% CREATE a rigidbody model and xyz object
rxyz = xyz.copy;
rxyz.model = rigidBody;
rxyz.data = xyz(:,rigidBody.ml,:);
rxyz.data = cat(1,permute(xyzh,[3,1,2]),rxyz.data);


mc = [0,0,1;...
      0,1,0;...
      0.5,0.2,0.5;...
      1,0,0];
ind = 1;
figure();hold('on');
for m = 1:4,
scatter3(rxyz(ind,m,1),rxyz(ind,m,2),rxyz(ind,m,3),30,mc(m,:),'filled')
end
daspect([1,1,1]);
grid('on');
% $$$ 
% $$$ 
% $$$ ang = create(MTADang,Trial,rxyz);


% CREATE an orthogonal coordinate system base on each Triad
markerTriad = compute_marker_triads(rxyz);

% COMPUTE the orthonormal triad basis for each frame
mtBasis = cf(@(x) bsxfun(@rdivide,x,sqrt(sum(x.^2,2))),...
             cf(@squeeze,...
                mat2cell(markerTriad.coor,...
                         ones([size(markerTriad.coor,1),1]),...
                         ones([size(markerTriad.coor,2),1]),3,3)...
               )...
);

goodIndex = 1;
reconstructedSolution ={};
for nck = 1:size(markerTriad.nck,1),
    goodTargetVector = sq([permute(neckxyz,[1,3,2])-rxyz(goodIndex,markerTriad.nck(nck,2),:)]);
    %goodTargetVector = sq([permute(nosexyz,[1,3,2])-rxyz(goodIndex,markerTriad.nck(nck,2),:)]);
    solutionBasisCoordinates = rref(cat(2,mtBasis{goodIndex,nck}',goodTargetVector));
    reconstructedSolution(nck,1) = {solutionBasisCoordinates(:,4)};
end


% $$$ ind =1000;
% $$$ figure,hold on
% $$$ plotSkeleton(Trial,xyz,ind);
% $$$ nck = 1;
% $$$ x = mtBasis{ind+1,nck}'*[reconstructedSolution{nck}]+sq(rxyz(ind+1,markerTriad.nck(nck,2),:));
% $$$ plot3(x(1),x(2),x(3),'*r')





% COMUPUTE the distance b
mrdist = nan([size(mtBasis,1),size(markerTriad.nck,1),1,2]);
sign = [1,-1];
markerTriad.solutions = nan([size(mtBasis,1),size(markerTriad.nck,1),size(rxyz,3)]);
markerTriad.oriMarkerNckCoor =nan([size(mtBasis,1),size(markerTriad.nck,1),1,size(rxyz,3)]);
% COMPUTE distance between reconstrution and original marker
for nck = 1:size(markerTriad.nck,1)
    markerTriad.solutions(:,nck,:) = ...
        cell2mat(cellfun(@transpose,...
                         cellfun(@(x,y) x'*[y], mtBasis(:,nck),... 
                                 repmat(reconstructedSolution(nck),...
                                        size(mtBasis,1),1),...
                                 'UniformOutput',false),...
                         'UniformOutput',false)...
                 );
end

%xyz       = preproc_xyz(Trial);
txyz = [];
for nck = 1:4
    txyz = cat(2,txyz,markerTriad.solutions(:,nck,:) + rxyz(:,markerTriad.nck(nck,2),:));
end
xyz.addMarker(['head_neck'],[200,255,200],{{'head_back','head_front',[0,0,1]}},mean(txyz(2:end,:,:),2));
%xyz.addMarker(['head_nose'],[.7,1,.7],{{'head_back','head_front',[0,0,1]}},mean(txyz(2:end,:,:),2));
xyz.addMarker(['head_nose_best'],[.7,1,.7],{{'head_back','head_front',[0,0,1]}},mean(txyz(2:end,:,:),2));


ind =10200;
figure,hold on
plotSkeleton(Trial,xyz,ind);
% $$$ nck = 1;
% $$$ x = mtBasis{ind+1,nck}'*reconstructedSolution{nck}+sq(rxyz(ind+1,markerTriad.nck(nck,2),:));
% $$$ plot3(x(1),x(2),x(3),'*r')


xyztrb = Trial.load('xyz','trb')
hcom = xyztrb.com(rigidBody);


xyz.addMarker(['head_trbcom'],[.7,1,.7],{{'head_back','head_front',[0,0,1]}},hcom);

ang = create(MTADang,Trial,xyz);

figure,histogram(ang(nniz(xyz),'head_neck','head_trbcom',3),linspace(0,30,100))

ind =10200;
figure,hold on
plotSkeleton(Trial,xyz,ind);


load(fullfile(Trial.spath,[Trial.filebase '.xyz-shift_fine.mat']));


ijk = cell(1,3);
[ijk{:}] = meshgrid(nj,ni,nk);

clist = 'rbgymck';


sp = [];
figure();hold('on');
mset = [1:7];
for m = mset
    %    sp(end+1) = subplot(2,4,m);
    svxyz = log10(sq(nvxyz(m,:,:,:)));
    iopts = [ijk,{svxyz},{prctile(svxyz(nniz(svxyz(:))),2)}];
    isos(m) = isosurface(iopts{:});
    p = patch(isos(m));
    iopts(end) = {p};
    isonormals(iopts{:});
    p.FaceColor = clist(m);
    p.EdgeColor = 'none';
    p.FaceAlpha = 0.3;
end
daspect([1 1 1]); 
camlight; lighting phong


% find the center of each field
mpos = zeros([7,3]);
mind = zeros([7,3]);
for m = 1:7,
    [mind(m,:),mv] = LocalMinimaN(sq(nvxyz(m,:,:,:)),100,100);
    mpos(m,:) = [nj(mind(m,2)),ni(mind(m,1)),nk(mind(m,3))];
end

nijk = cf(@(x)  x(:),  ijk);
nijk = cat(2,nijk{:});

varLines = [];
for m = mset
% $$$     distToMv = sqrt(sum(bsxfun(@minus,nijk,mpos(m,:)).^2,2));
% $$$     tnvxyz = sq(nvxyz(m,:,:,:));
% $$$     tnvxyz = 1./(tnvxyz(:));
% $$$     %dind = distToMv<30;
% $$$     dind = tnvxyz(:)< prctile(tnvxyz(nniz(tnvxyz(:))),20);
% $$$     [~,S,V] = svd(bsxfun(@times,tnvxyz(dind,:),nijk(dind,:)),0);
    [~,S,V] = svd(bsxfun(@minus,isos(m).vertices,mpos(m,:)),0);
    for vind = 1
        %axes(sp(m));
        quiver3(mpos(m,1),mpos(m,2),mpos(m,3),V(1,vind),V(2,vind),V(3,vind),30)
    end
    varLines(m,:) = V(:,1);
end

scp = nan([7,7,3]);
tcq = nan([7,7,3]);
for p = 1:7,
    for q = p+1:7,
u = varLines(p,:);
v = varLines(q,:);
w = mpos(p,:)-mpos(q,:);
a = dot(u,u);
b = dot(u,v);
c = dot(v,v);
d = dot(u,w);
e = dot(v,w);

sc = (b*e-c*d)/(a*c-b^2);
tc = (a*e-b*d)/(a*c-b^2);

scp(p,q,:) = mpos(p,:)+sc*varLines(p,:);
tcq(p,q,:) = mpos(q,:)+tc*varLines(q,:);
end
end
% $$$ plot3(scp(1),scp(2),scp(3),'*b')
% $$$ plot3(tcq(1),tcq(2),tcq(3),'*m')


mpoint = nanmean(nanmean(cat(3,reshape(scp,[],3),reshape(tcq,[],3)),3));

scatter3(mpoint(1),mpoint(2),mpoint(3),50,[0,0,1],'filled')



xyz.addMarker(['head_neck_est'],[255,0,0],{{'head_back','head_front',[1,0,1]}},bsxfun(@plus,nx*mpoint(2)+ny*mpoint(1)+nz*mpoint(3),xyz(:,'hcom',:)));

i = 1;
xyz.addMarker(['head_neck_est'],[20,20,255],{{'head_back','head_front',[1,0,1]}},bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'hcom',:)));


ind = 20000;
figure,hold on
plotSkeleton(Trial,xyz,ind);
scatter3(xyz(ind,'hcom',1),xyz(ind,'hcom',2),xyz(ind,'hcom',3),30,[1,0,1],'filled');
scatter3(xyz(ind,'head_center',1),xyz(ind,'head_center',2),xyz(ind,'head_center',3),30,[1,0,1],'filled');
scatter3(xyz(ind,'head_centerS',1),xyz(ind,'head_centerS',2),xyz(ind,'head_centerS',3),30,[1,0,1],'filled');





xyz.data(:,xyz.model.gmi('head_back'),:) =  bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_back',:));
xyz.data(:,xyz.model.gmi('head_left'),:) =  bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_left',:));
xyz.data(:,xyz.model.gmi('head_front'),:) =  bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_front',:));
xyz.data(:,xyz.model.gmi('head_right'),:) =  bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'head_right',:));


