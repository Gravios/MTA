OwnDir = '/storage/gravio/nextcloud/';
FigDir = 'Shared/Behavior Paper/Figures/Suplementary';
mkdir(fullfile(OwnDir,FigDir));



Trial = MTATrial.validate('Ed03-20140624.cof.all');


% stereotaxic coordinates of head marker positions
xyzh = [ 23.71, 14.43, 74.78;... head_back
         52.44,-11.00, 66.56;... head_left
         76.82,  4.47, 65.67;... head_front
         51.46, 35.80, 65.67;... head_right
];

% $$$ nosexyz = [79.47, 10.00, 31.12];
% $$$ neckxyz = [79.47-30, 10.00, 31.12];

yshift  = min(xyzh(:,2));
xyzh = bsxfun(@times,...
              bsxfun(@plus,xyzh,[0,abs(yshift),0]),...
              [1,-1,1]);

nosexyz = [79.47, 10.00, 31.12];
neckxyz = [79.47-40, 10.00, 31.12];


newPoints = {};
newPoints{1} = ([79.47, 10.00, 31.12]+[0,abs(yshift),0]).*[1,-1,1];
newPoints{2} = ([79.47-40, 10.00, 31.12]+[0,abs(yshift),0]).*[1,-1,1];
newPointsLabels = {'head_nose','head_neck'};


% Ed05
Trial = MTATrial.validate('Ed05-20140529.ont.all');
xyzh = [ 22.46, 11.96, 61.07;... head_back
         52.42,-18.00, 60.24;... head_left                  
         80.00, -2.38, 57.25;... head_front
         51.48, 34.50, 65.02;... head_right
];

nosexyz = [85.68, 9.65, 32.39];
neckxyz = [85.68-30, 9.65, 32.39];


% Correct with rotation?
% $$$ m = 3;
% $$$ ax_ord = [1,2,3];
% $$$ j =1:3;
% $$$ head_norm = bsxfun(@rdivide,sq(evec(ind,m,ax_ord)),sqrt(sum(evec(ind,m,ax_ord).^2,3)));
% $$$ head_kron = reshape(repmat(head_norm',3,1).*head_norm(:,j(ones(3,1),:)).',[3,3,size(head_norm,1)]);
% $$$ j = [ 0,-1, 1;...
% $$$       1, 0,-1;...
% $$$       -1, 1, 0];
% $$$ k = [1,3,2;...
% $$$      3,1,1;...
% $$$      2,1,1];
% $$$ head_cpm = reshape(head_norm(:,k)',3,3,size(head_norm,1)).*repmat(j,[1,1,size(head_norm,1)]);
% $$$ 
% $$$ j =1:3;  
% $$$ rot_ang = deg2rad(rots(d));
% $$$ head_rotMat = cos(rot_ang)*repmat(eye(3),[1,1,size(head_norm,1)])...
% $$$     +sin(rot_ang)*head_cpm...
% $$$     +(1-cos(rot_ang))*head_kron;


yshift  = min(xyzh(:,2));
xyzh = bsxfun(@times,...
              bsxfun(@plus,xyzh,[0,abs(yshift),0]),...
              [1,-1,1]);

newPoints = {};
newPoints{1} = ([85.68   , 9.65, 32.39]+[0,abs(yshift),0]).*[1,-1,1];
newPoints{2} = ([85.68-40, 9.65, 32.39]+[0,abs(yshift),0]).*[1,-1,1];
newPointsLabels = {'head_nose','head_neck'};



% LOAD xyz data
xyz       = preproc_xyz(Trial);


nz = -cross(xyz(:,'head_back',:)-xyz(:,'hcom',:),xyz(:,'head_left',:)-xyz(:,'hcom',:));
nz = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3)));
nm = nz.*20+xyz(:,'hcom',:);
xyz.addMarker('htx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);

ny = cross(xyz(:,'htx',:)-xyz(:,'hcom',:),xyz(:,'head_back',:)-xyz(:,'hcom',:));
ny = bsxfun(@rdivide,ny,sqrt(sum((ny).^2,3)));
nm = ny.*20+xyz(:,'hcom',:);
xyz.addMarker('hrx',[128,255,128],{{'head_back','head_front',[0,0,1]}},nm);

nx = cross(xyz(:,'hrx',:)-xyz(:,'hcom',:),xyz(:,'htx',:)-xyz(:,'hcom',:));
nx = bsxfun(@rdivide,nx,sqrt(sum((nx).^2,3)));
nm = nx.*20+xyz(:,'hcom',:);    
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
scatter3(newPoints{1}(1),newPoints{1}(2),newPoints{1}(3),30,mc(m,:),'filled')
daspect([1,1,1]);
grid('on');


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


for p = 1:numel(newPoints)
    goodIndex = 1;
    reconstructedSolution ={};
    for nck = 1:size(markerTriad.nck,1),
        goodTargetVector = sq([permute(newPoints{p},[1,3,2])-rxyz(goodIndex,markerTriad.nck(nck,2),:)]);
        solutionBasisCoordinates = rref(cat(2,mtBasis{goodIndex,nck}',goodTargetVector));
        reconstructedSolution(nck,1) = {solutionBasisCoordinates(:,4)};
    end


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

    txyz = [];
    for nck = 1:4
        txyz = cat(2,txyz,markerTriad.solutions(:,nck,:) + rxyz(:,markerTriad.nck(nck,2),:));
    end
    xyz.addMarker(newPointsLabels{p},[200,50,100],{{'head_back','head_front',[0,0,1]}},mean(txyz(2:end,:,:),2));

end

xyzsehs = Trial.load('xyz','sehs');

ind =1200;
figure,hold on
subplot2(4,8,[1:3],[1:4]);
plotSkeleton(Trial,xyz,ind);
plotSkeleton(Trial,xyzsehs,ind);
title('trb estimation evaluation against emperical measurements');
xyz.addMarker(['head_trbcom'],[.7,1,.7],{{'head_back','head_front',[0,0,1]}},xyzsehs(:,'hcom',:));
ang = create(MTADang,Trial,xyz);
subplot2(4,8,4,[2:3]);
histogram(ang(nniz(xyz),'head_neck','head_trbcom',3),linspace(0,20,100)),axis tight
ylabel('count')
xlabel('distance between neck position estimates (mm)')




% PLOT the info for trb estimation 
load(fullfile(Trial.spath,[Trial.filebase '.xyz-shift_fine.mat']));


ijk = cell(1,3);
[ijk{:}] = meshgrid(nj,ni,nk);

clist = 'rbgymck';



%figure();
subplot2(4,8,1:4,6:8);
hold('on');

mset = [1:7];
for m = mset
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

saveas(gcf,fullfile(OwnDir,FigDir,['trb_estimation_confirmation']),'fig')
print(gcf,'-depsc2',fullfile(OwnDir,FigDir,['trb_estimation_confirmation.eps']));
print(gcf,'-dpng',fullfile(OwnDir,FigDir,['trb_estimation_confirmation.png']));


xyz.addMarker(['head_neck_est'],[255,0,0],{{'head_back','head_front',[1,0,1]}},bsxfun(@plus,nx*mpoint(2)+ny*mpoint(1)+nz*mpoint(3),xyz(:,'hcom',:)));

i = 1;
xyz.addMarker(['head_neck_est'],[20,20,255],{{'head_back','head_front',[1,0,1]}},bsxfun(@plus,nx*ni(mind(i,1))+ny*nj(mind(i,2))+nz*nk(mind(i,3)),xyz(:,'hcom',:)));







transform_rigidBody(Trial);
% diagnostics
nxyz.data(~nniz(xyz),:,:) = eps;
nxyz.data(~nniz(xyz),:,:) = 0;
nxyz.save;

Trial = MTATrial.validate('Ed05-20140529.ont.all');
Trial = MTATrial.validate('Ed05-20140529.ont.all');
xyztrb = Trial.load('xyz','trb');
xyz = Trial.load('xyz');

ind = 20000;
figure,hold on
plotSkeleton(Trial,xyztrb,ind);
plotSkeleton(Trial,xyz,ind);

angtrb = create(MTADang,Trial,xyztrb);
figure, plot(nunity(angtrb(:,'spine_upper','head_neck',3)))
hold on,plot(nunity(angtrb(:,'spine_upper','head_neck',2)))



Trials = af(@(t)  MTATrial.validate(t),                        get_session_list('hand_labeled'));
t = 6;
transform_rigidBody(Trials{t},true,false);
fet_spline_spine(Trials{t},'overwrite',true)
xyz = preproc_xyz(Trials{t},'SPLINE_SPINE_HEAD_EQI',true);

load_normalization_parameters_mapminmax('fet_mis','overwrite',true);


Trials = af(@(t)  MTATrial.validate(t),  get_session_list('BHV_S4H5'));
Trials([1,5,7,8,10,27]) = [];



t = 21;
transform_rigidBody(Trials{t},true,true);  
fet_spline_spine(Trials{t},'overwrite',true)
preproc_xyz(Trials{t},'SPLINE_SPINE_HEAD_EQI',true)

cf(@(t)  transform_rigidBody(t,true,true),           Trials);
cf(@(t)  fet_spline_spine(t,'overwrite',true),        Trials);
cf(@(t)  preproc_xyz(t,'SPLINE_SPINE_HEAD_EQI',true), Trials);