

Trial = MTATrial.validate('jg05-20120317.cof.all');

xyz = infer_virtual_joint_from_rigidbody_kinematics(Trial);

xyz = Trial.load('xyz','crb');

%% -- csv time ---
mxyz = importdata('/storage/antsiro/data/lab/fabian/dbase/FS11_20220322-185620/FS11_20220322-185620_actual.csv');

% 4 Markers 4 dims
43:4:58

nxyz = zeros([size(mxyz.data,1),4,3]);
for d = 1:3,
    nxyz(:,:,d) = mxyz.data(:,42+d:4:58);
end

nxyz(:,:,[2,3]) = nxyz(:,:,[3,2]);

nxyz = nxyz*1000;

fxyz = nxyz;
fxyz(nniz(nxyz),:,:) = ButFilter(nxyz(nniz(nxyz),:,:),4,20/(120/2),'low');

hcom = mean(fxyz,2);

nz = -cross(fxyz(:,1,:)-hcom,fxyz(:,2,:)-hcom);
nz = bsxfun(@rdivide,nz,sqrt(sum((nz).^2,3))); 
ny = cross(nz,fxyz(:,1,:)-hcom);
ny = bsxfun(@rdivide,ny,sqrt(sum((ny).^2,3)));
nx = cross(nz,ny);
nx = bsxfun(@rdivide,nx,sqrt(sum((nx).^2,3)));
rigidBodyBasis = cat(3,sq(nx),sq(ny),sq(nz));
headCenterOfMass = sq(hcom);

sampleRate = 120;

headSpeed = sqrt(sum(clip(multiprod(rigidBodyBasis,circshift(sq(headCenterOfMass),-1)-circshift(sq(headCenterOfMass),1),[2,3],[2]).*sampleRate/10/2,-1000,1000).^2,2,'omitnan'));



% GET "appropriate" subset of of data for fitting
nind = ( headSpeed  > 5 )  &  ( headSpeed  < 80 );%  &  ( headCenterOfMass(:,3) < 150+450 ); 
% Units              cm/s                    cm/s                                 mm

shifts = [-100:100]; 

virtualJointCenter = zeros([1,size(xyz,3)]);
for dim = 1:size(xyz,3)
    offset = [0,0,0];
    offset(dim) = 1;
    grad = [];         
    for shift = shifts,
        scom = headCenterOfMass + multiprod(rigidBodyBasis,offset*shift,[2,3],2);
        dcom = circshift(scom,-1)-circshift(scom,1);
        scomVelocity = sum(clip(multiprod(rigidBodyBasis(nind,:,:),dcom(nind,:),[2,3],[2]),-100,100).^2,'omitnan');
        grad = cat(1,grad,log10(scomVelocity));
    end
    [~,mind] = min(grad);
    virtualJointCenter(1,dim) = mean(shifts(mind));
    figure,plot(shifts,grad);title(['dim ' num2str(dim)]);
end




figure();
hold('on');
mclr = 'brgm';
ind = 4000:4200;
for m = 1:4
plot3(fxyz(ind,m,1),fxyz(ind,m,2),fxyz(ind,m,3),['-',mclr(m)]);
end