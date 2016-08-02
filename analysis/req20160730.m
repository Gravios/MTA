


dir = '/storage/nickdg/data/VR Engagement/';
data = load(fullfile(dir,'justin_VRObj_ratdata.mat'));
ndata = data.data;
%data.data(:,1): time

% rigidbody center of mass
%data.data(:,2): Y    Nick's X
%data.data(:,3): Z    Nick's Y
%data.data(:,4): X    Nick's Z

% euler angles of rigid body relative to room cooridate system
%data.data(:,5): Roll
%data.data(:,6): Pitch
%data.data(:,7): Yaw

% Some orientation vector nick made
%data.data(:,8): Y
%data.data(:,9): Z
%data.data(:,10): X

% Session Id
%data.data(:,11): session id


% euler ang vector
eang = deg2rad(data.data(:,[5,6,7]));

yaw = eang(:,1);
pitch = eang(:,2);
roll = eang(:,3);

% $$$ rMat = cell([size(eang,1),1]);
% $$$ for i = 1:size(eang,1),
% $$$ rMat{i} = [cos(yaw(i)).*cos(pitch(i)),...
% $$$        -cos(yaw(i)).*sin(pitch(i)).*sin(roll(i))-sin(yaw(i)).*cos(roll(i)),...
% $$$        -cos(yaw(i)).*sin(pitch(i)).*cos(roll(i))+sin(yaw(i)).*sin(roll(i));...
% $$$         ...
% $$$         sin(yaw(i)).*cos(pitch(i)),...
% $$$        -sin(yaw(i)).*sin(pitch(i)).*sin(roll(i))+cos(yaw(i)).*cos(roll(i)),...
% $$$        -sin(yaw(i)).*sin(pitch(i)).*cos(roll(i))-cos(yaw(i)).*sin(roll(i));...
% $$$         ...
% $$$         sin(pitch(i)),...
% $$$         cos(pitch(i)).*sin(roll(i)),...
% $$$         cos(pitch(i)).*sin(roll(i))];
% $$$ end
% $$$ 
% $$$ rMat = cell([size(eang,1),1]);
% $$$ for i = 1:size(eang,1),
% $$$ rMat{i} = [cos(roll(i)).*cos(pitch(i)),...
% $$$            cos(roll(i)).*sin(pitch(i)).*sin(yaw(i))-sin(roll(i)).*cos(yaw(i)),...
% $$$            cos(roll(i)).*sin(pitch(i)).*cos(yaw(i))+sin(roll(i)).*sin(yaw(i));...
% $$$                ...
% $$$            sin(roll(i)).*cos(pitch(i)), ...
% $$$            sin(roll(i)).*sin(pitch(i)).*sin(yaw(i))+cos(roll(i)).*cos(yaw(i)),...
% $$$            sin(roll(i)).*sin(pitch(i)).*cos(yaw(i))-cos(roll(i)).*sin(yaw(i));...
% $$$                ...
% $$$           -sin(roll(i)),...
% $$$            cos(pitch(i)).*sin(yaw(i)),...
% $$$            cos(pitch(i)).*cos(yaw(i))];
% $$$ end

% RxRyRz
rMat = cell([size(eang,1),1]);
for i = 1:size(eang,1),
rMat{i} = [cos(pitch(i)).*cos(yaw(i)),...
          -cos(pitch(i)).*sin(yaw(i)),...
           sin(pitch(i));...
               ...
           cos(roll(i)).*sin(yaw(i))+sin(roll(i)).*sin(pitch(i)).*cos(yaw(i)),...
           cos(roll(i)).*cos(yaw(i))-sin(roll(i)).*sin(pitch(i)).*sin(yaw(i)),...
          -sin(roll(i)).*cos(pitch(i));...
               ...
          sin(roll(i)).*sin(yaw(i))-cos(roll(i)).*sin(pitch(i)).*cos(yaw(i)),...
          sin(roll(i)).*cos(yaw(i))+cos(roll(i)).*sin(pitch(i)).*sin(yaw(i)),...
          cos(roll(i)).*cos(pitch(i))];
end



blen = 0.05;
vt = cellfun(@mtimes,rMat,repmat({[blen;0;0]},[size(eang,1),1]),'UniformOutput',false);
evec(:,1,:) = reshape(cell2mat(vt'),3,[])';
vt = cellfun(@mtimes,rMat,repmat({[0;blen;0]},[size(eang,1),1]),'UniformOutput',false);
evec(:,2,:) = reshape(cell2mat(vt'),3,[])';
vt = cellfun(@mtimes,rMat,repmat({[0;0;blen]},[size(eang,1),1]),'UniformOutput',false);
evec(:,3,:) = reshape(cell2mat(vt'),3,[])';

evec = evec(:,:,[3,1,2]);


% Get the session ids (30 total)
sessionIds = unique(data.data(:,11));

xyzSampleRate = 150;
% store XYZ rigid body center of mass in meters
% X -> long axis of maze
% Y -> short axis of maze
% Z -> Height
% xyz(time,marker,dimension)
xyz(:,1,:) = data.data(:,[4,2,3]);
% Add a low passed (0.8 Hz) version of center of mass with slight time shift
xyz(:,2,:) = circshift(ButFilter(xyz(:,1,:),3,[.8]./(xyzSampleRate/2),'low'),round(xyzSampleRate/4));
% Add a low passed (0.2 Hz) version of center of mass with large time shift
xyz(:,3,:) = circshift(ButFilter(xyz(:,1,:),3,[.2]./(0.5*xyzSampleRate),'low'),round(xyzSampleRate));
% Add orientation vector + center of mass reduce in length
xyz(:,4,:) = 0.05.*data.data(:,[10,8,9])+data.data(:,[4,2,3]);
% Add a low passed (2 Hz) version of center of mass with no time shift
xyz(:,5,:) = ButFilter(xyz(:,1,:),3,[2]./(xyzSampleRate/2),'low');
% Add markers surrounding the center of mass to create a head basis
xyz(:,6,:) =  sq(evec(:,1,:))+data.data(:,[4,2,3]);
xyz(:,7,:) =  sq(evec(:,2,:))+data.data(:,[4,2,3]);
xyz(:,8,:) =  sq(evec(:,3,:))+data.data(:,[4,2,3]);


figure,hold on
ind = 1000;
plot3(xyz(ind,1,1),xyz(ind,1,2),xyz(ind,1,3),'.m')
plot3(xyz(ind,6,1),xyz(ind,6,2),xyz(ind,6,3),'.r')
plot3(xyz(ind,4,1),xyz(ind,4,2),xyz(ind,4,3),'.c')
plot3(xyz(ind,7,1),xyz(ind,7,2),xyz(ind,7,3),'.g')
plot3(xyz(ind,8,1),xyz(ind,8,2),xyz(ind,8,3),'.k')
daspect([1,1,1])

% create low passed (2.4 Hz) speed in xy plane
markers = [1,5];
vxy = sqrt(sum([zeros([1,numel(markers),2]);diff(xyz(:,markers,[1,2]))].^2,3)).*xyzSampleRate;



% markerDiffMatrix
% diffMat(time,marker1,marker2,dim) 
% create matrix where each marker's position is subtracted from
% every other marker

nframe = size(xyz,1); %nframe: number of frames (time)
nmar   = size(xyz,2);  %nmar: number markers 
ndim   = size(xyz,3); %ndim: number of spatial dimensions (xyz)
j =1:nmar;
diffMat = permute(cat(4,permute(reshape(repmat(xyz(:,:,1)',nmar,1)-xyz(:,j(ones(nmar,1),:),1).',[nmar,nmar,nframe]),[3,1,2]),permute(reshape(repmat(xyz(:,:,2)',nmar,1)-xyz(:,j(ones(nmar,1),:),2).',[nmar,nmar,nframe]),[3,1,2]),permute(reshape(repmat(xyz(:,:,3)',nmar,1)-xyz(:,j(ones(nmar,1),:),3).',[nmar,nmar,nframe]),[3,1,2])),[1,3,2,4]);


% Create Angles between all markers
% ang(time,marker1,marker2,dim) 
%     dim: spherical coordinates (theta, phi   , rho     )
%                                (yaw  , pitch , distance)
ang = zeros(size(xyz,1),size(xyz,2),size(xyz,2),3);
for i=1:size(xyz,2),
    for j=1:size(xyz,2),
        if i==j,continue,end                    
        switch size(xyz,3)
          case 3
            tang =cell(1,3);
            [tang{:}] = cart2sph(diffMat(:,i,j,1),diffMat(:,i,j,2),diffMat(:,i,j,3));
          case 2
            tang =cell(1,2);
            [tang{:}] = cart2pol(diffMat(:,i,j,1),diffMat(:,i,j,2),diffMat(:,i,j,3));
        end
        ang(:,i,j,:) = cell2mat(tang);
    end
end
ang(ang(:,1,2,2)~=0,1,1,1)=1;



%% Find periods of walking with just the head
%
% Come back to this later maybe
%
% $$$ cang = circ_dist(ang(:,2,3,1),ang(:,1,3,1));
% $$$ dcang = circ_dist(circshift(cang,10),circshift(cang,-10));
% $$$ nn
% $$$ bang = clip(ang(:,2,3,3)./sqrt(sq(sum(GetSegs(circshift(dcang,25),1:size(dcang),50).^2)))'./10,0,1000);
% $$$ bang(isnan(bang)|bang==0|isinf(bang)) = 1e-10;
% $$$ bang = ButFilter(bang,3,2.4./(0.5*xyzSampleRate),'low');
% $$$ 
% $$$ bThresh = 0.01;
% $$$ mDur = xyzSampleRate*2;
% $$$ 
% $$$ wper = ThreshCross(bang,bThresh,mDur);
% $$$ 
% $$$ 
% $$$ w = 1;
% $$$ ind = wper(w,1):wper(w,2);
% $$$ 
% $$$ figure,plot3(xyz(ind,1,1),xyz(ind,1,2),xyz(ind,1,3))
% $$$ figure,plot(ang(ind,1,4,1))
% $$$ 
% $$$ angMean = ones([size(wper,1),1]);
% $$$ angVar  = ones([size(wper,1),1]);
% $$$ 
% $$$ for i = 1:size(wper,1),
% $$$     ind = wper(i,1)+50:wper(i,2)-50;    
% $$$     tmp = bsxfun(@circ_dist,ang(ind,1,4,1)',circshift(ang(ind,2,3,1),-round(0.75*xyzSampleRate)));
% $$$     angMean(i) = circ_mean(tmp(:));
% $$$     angVar(i)  = circ_var(tmp(:));  
% $$$ end
% $$$ 
% $$$ figure,plot3(xyz(ind,1,1),xyz(ind,1,2),xyz(ind,1,3))
% $$$ hold on,plot3(xyz(ind,2,1),xyz(ind,2,2),xyz(ind,2,3))
% $$$ 
% $$$ figure,hist(angMean,100)
% $$$ Lines(circ_median(angMean),[],'r')



% Sniffing
figure,plot([1:size(ang,1)]./xyzSampleRate,[0;diff(ButFilter(ang(:,6,5,3),3,20/(0.5.*xyzSampleRate),'low'))])

% Sniffing
figure,plot([1:size(ang,1)]./xyzSampleRate,[0;diff(ButFilter(ang(:,4,5,3),3,20/(0.5.*xyzSampleRate),'low'))])


% Find Best orientation (good luck)


ind = data.data(:,11)==sessionIds(2);

figure,hold on
moffset = 1200;
sind = find(ind,1,'first')+moffset;
plot3(xyz(sind,1,1),xyz(sind,1,2),xyz(sind,1,3),'.m')
plot3(xyz(sind,6,1),xyz(sind,6,2),xyz(sind,6,3),'.r')
plot3(xyz(sind,7,1),xyz(sind,7,2),xyz(sind,7,3),'.g')
plot3(xyz(sind,8,1),xyz(sind,8,2),xyz(sind,8,3),'.k')


spw = [];
tots = 0:5:380;
rots = 5:5:180;
for d =  1:numel(rots)
    m = 3;
    ax_ord = [1,2,3];
    j =1:3;
    head_norm = bsxfun(@rdivide,sq(evec(ind,m,ax_ord)),sqrt(sum(evec(ind,m,ax_ord).^2,3)));
    head_kron = reshape(repmat(head_norm',3,1).*head_norm(:,j(ones(3,1),:)).',[3,3,size(head_norm,1)]);
    j = [ 0,-1, 1;...
          1, 0,-1;...
          -1, 1, 0];
    k = [1,3,2;...
         3,1,1;...
         2,1,1];
    head_cpm = reshape(head_norm(:,k)',3,3,size(head_norm,1)).*repmat(j,[1,1,size(head_norm,1)]);

    j =1:3;  
    rot_ang = deg2rad(rots(d));
    head_rotMat = cos(rot_ang)*repmat(eye(3),[1,1,size(head_norm,1)])...
        +sin(rot_ang)*head_cpm...
        +(1-cos(rot_ang))*head_kron;

    % Rotated marker;
    ovec = evec(ind,1,:);

    nvec = permute(sum(head_rotMat.*permute(reshape(ovec(:,j(ones(3,1),:)),[size(head_norm,1),3,3]),[2,3,1]),2),[3,2,1]);

    nxyz = cat(2,xyz(ind,[6:8],:),nvec);

    nvec = cat(2,evec(ind,:,:),nvec);

    m = 4;
    ax_ord = [1,2,3];
    j =1:3;
    head_norm = bsxfun(@rdivide,sq(nvec(:,m,ax_ord)),sqrt(sum(nvec(:,m,ax_ord).^2,3)));
    head_kron = reshape(repmat(head_norm',3,1).*head_norm(:,j(ones(3,1),:)).',[3,3,size(head_norm,1)]);
    j = [ 0,-1, 1;...
          1, 0,-1;...
          -1, 1, 0];
    k = [1,3,2;...
         3,1,1;...
         2,1,1];
    head_cpm = reshape(head_norm(:,k)',3,3,size(head_norm,1)).*repmat(j,[1,1,size(head_norm,1)]);


    for t =  1:numel(tots)
        j =1:3;  
        rot_ang = deg2rad(tots(t));
        head_rotMat = cos(rot_ang)*repmat(eye(3),[1,1,size(head_norm,1)])...
            +sin(rot_ang)*head_cpm...
            +(1-cos(rot_ang))*head_kron;

        % Rotated marker;
        %ovec = nvec(ind,2,:);        
        ovec = nvec(:,2,:);
        ovec = cross(nvec(:,2,:),nvec(:,4,:));
        ovec = bsxfun(@rdivide,ovec,sqrt(sum(ovec.^2,3))).*blen;


        tmark = permute(sum(head_rotMat.*permute(reshape(ovec(:,j(ones(3,1),:)),[size(head_norm,1),3,3]),[2,3,1]),2),[3,1,2])+data.data(ind,[4,2,3]);

        sxyz = cat(2,xyz(ind,5,:),permute(tmark,[1,3,2]));
        

        nframe = size(sxyz,1); %nframe: number of frames (time)
        nmar   = size(sxyz,2);  %nmar: number markers 
        ndim   = size(sxyz,3); %ndim: number of spatial dimensions (sxyz)
        j =1:nmar;
        diffMat = permute(cat(4,permute(reshape(repmat(sxyz(:,:,1)',nmar,1)-sxyz(:,j(ones(nmar,1),:),1).',[nmar,nmar,nframe]),[3,1,2]),permute(reshape(repmat(sxyz(:,:,2)',nmar,1)-sxyz(:,j(ones(nmar,1),:),2).',[nmar,nmar,nframe]),[3,1,2]),permute(reshape(repmat(sxyz(:,:,3)',nmar,1)-sxyz(:,j(ones(nmar,1),:),3).',[nmar,nmar,nframe]),[3,1,2])),[1,3,2,4]);

        ang = zeros(size(sxyz,1),size(sxyz,2),size(sxyz,2),3);
        for i=1:size(sxyz,2),
            for j=1:size(sxyz,2),
                if i==j,continue,end                    
                switch size(sxyz,3)
                  case 3
                    tang =cell(1,3);
                    [tang{:}] = cart2sph(diffMat(:,i,j,1),diffMat(:,i,j,2),diffMat(:,i,j,3));
                  case 2
                    tang =cell(1,2);
                    [tang{:}] = cart2pol(diffMat(:,i,j,1),diffMat(:,i,j,2),diffMat(:,i,j,3));
                end
                ang(:,i,j,:) = cell2mat(tang);
            end
        end
        ang(ang(:,1,2,2)~=0,1,1,1)=1;
        
        bang = ButFilter(diff(ButFilter(ang(:,1,2,3),3,20./(0.5*xyzSampleRate),'low')),3,[6,15]./(0.5*xyzSampleRate),'bandpass');
        spw(d,t) = nanmedian(sqrt(sum(GetSegs(bang,1:round(xyzSampleRate/2):size(bang,1)-round(xyzSampleRate),round(xyzSampleRate)).^2)));

        plot3(tmark(moffset,1),tmark(moffset,2),tmark(moffset,3),'*r')
    end
end


% spw is the map where the peaks represents the front and back of
% the head
% the indicies (d,t) are the two angle thingies you need to get the final
% orientation vector. 
%
% To Do: automatic peak detection and functionalize the code above
%        
% For now: find the two indices plug them in above and get the new
%          marker called tmark as the new orientation vector.

figure,imagesc(spw')

nspw  = spw(nniz(spw),:);
figure,imagesc(nspw')

[mins,minv] = LocalMinimaN(-nspw,-mean(nspw(:)),20);

nrots = rots(nniz(spw));
ntots = tots;

% Get best vector set
bvecs = [];

for dt = mins'
    m = 3;
    ax_ord = [1,2,3];
    j =1:3;
    head_norm = bsxfun(@rdivide,sq(evec(ind,m,ax_ord)),sqrt(sum(evec(ind,m,ax_ord).^2,3)));
    head_kron = reshape(repmat(head_norm',3,1).*head_norm(:,j(ones(3,1),:)).',[3,3,size(head_norm,1)]);
    j = [ 0,-1, 1;...
          1, 0,-1;...
          -1, 1, 0];
    k = [1,3,2;...
         3,1,1;...
         2,1,1];
    head_cpm = reshape(head_norm(:,k)',3,3,size(head_norm,1)).*repmat(j,[1,1,size(head_norm,1)]);


    j =1:3;  
    rot_ang = deg2rad(nrots(dt(1)));
    head_rotMat = cos(rot_ang)*repmat(eye(3),[1,1,size(head_norm,1)])...
        +sin(rot_ang)*head_cpm...
        +(1-cos(rot_ang))*head_kron;

    % Rotated marker;
    ovec = evec(ind,1,:);

    nvec = permute(sum(head_rotMat.*permute(reshape(ovec(:,j(ones(3,1),:)),[size(head_norm,1),3,3]),[2,3,1]),2),[3,2,1]);

    nxyz = cat(2,xyz(ind,[6:8],:),nvec);

    nvec = cat(2,evec(ind,:,:),nvec);

    m = 4;
    ax_ord = [1,2,3];
    j =1:3;
    head_norm = bsxfun(@rdivide,sq(nvec(:,m,ax_ord)),sqrt(sum(nvec(:,m,ax_ord).^2,3)));
    head_kron = reshape(repmat(head_norm',3,1).*head_norm(:,j(ones(3,1),:)).',[3,3,size(head_norm,1)]);
    j = [ 0,-1, 1;...
          1, 0,-1;...
          -1, 1, 0];
    k = [1,3,2;...
         3,1,1;...
         2,1,1];
    head_cpm = reshape(head_norm(:,k)',3,3,size(head_norm,1)).*repmat(j,[1,1,size(head_norm,1)]);


    j =1:3;  
    rot_ang = deg2rad(ntots(dt(2)));
    head_rotMat = cos(rot_ang)*repmat(eye(3),[1,1,size(head_norm,1)])...
        +sin(rot_ang)*head_cpm...
        +(1-cos(rot_ang))*head_kron;

    % Rotated marker;
    %ovec = nvec(ind,2,:);        
    ovec = nvec(:,2,:);
    ovec = cross(nvec(:,2,:),nvec(:,4,:));
    ovec = bsxfun(@rdivide,ovec,sqrt(sum(ovec.^2,3))).*blen;

    tmark = permute(permute(sum(head_rotMat.*permute(reshape(ovec(:,j(ones(3,1),:)),[size(head_norm,1),3,3]),[2,3,1]),2),[3,1,2])+data.data(ind,[4,2,3]),[1,3,2]);

    bvecs = cat(2,bvecs,tmark);
end


cxyz = cat(2,xyz(ind,[1,2,4],:),bvecs,xyz(ind,[6,7,8],:));

% markerDiffMatrix
% diffMat(time,marker1,marker2,dim) 
% create matrix where each marker's position is subtracted from
% every other marker

nframe = size(cxyz,1); %nframe: number of frames (time)
nmar   = size(cxyz,2);  %nmar: number markers 
ndim   = size(cxyz,3); %ndim: number of spatial dimensions (xyz)
j =1:nmar;
diffMat = permute(cat(4,permute(reshape(repmat(cxyz(:,:,1)',nmar,1)-cxyz(:,j(ones(nmar,1),:),1).',[nmar,nmar,nframe]),[3,1,2]),permute(reshape(repmat(cxyz(:,:,2)',nmar,1)-cxyz(:,j(ones(nmar,1),:),2).',[nmar,nmar,nframe]),[3,1,2]),permute(reshape(repmat(cxyz(:,:,3)',nmar,1)-cxyz(:,j(ones(nmar,1),:),3).',[nmar,nmar,nframe]),[3,1,2])),[1,3,2,4]);


% Create Angles between all markers
% ang(time,marker1,marker2,dim) 
%     dim: spherical coordinates (theta, phi   , rho     )
%                                (yaw  , pitch , distance)
ang = zeros(size(cxyz,1),size(cxyz,2),size(cxyz,2),3);
for i=1:size(cxyz,2),
    for j=1:size(cxyz,2),
        if i==j,continue,end                    
        switch size(cxyz,3)
          case 3
            tang =cell(1,3);
            [tang{:}] = cart2sph(diffMat(:,i,j,1),diffMat(:,i,j,2),diffMat(:,i,j,3));
          case 2
            tang =cell(1,2);
            [tang{:}] = cart2pol(diffMat(:,i,j,1),diffMat(:,i,j,2),diffMat(:,i,j,3));
        end
        ang(:,i,j,:) = cell2mat(tang);
    end
end
ang(ang(:,1,2,2)~=0,1,1,1)=1;

mang = [];
for m = 1:numel(minv),
    vind = vxy(ind,2)>0.2;
    mang(m) = circ_mean(circ_dist(ang(vind,2,1,1),ang(vind,1,m+3,1)));
    figure,rose(circ_dist(ang(vind,2,1,1),ang(vind,1,m+3,1)),100)
end


ndata(ind,[10,8,9]) = bvecs(:,min(abs(mang))==abs(mang),:);


% $$$ 
% $$$ figure,hold on
% $$$ sind = 1200;
% $$$ plot3(xyz(sind,1,1),xyz(sind,1,2),xyz(sind,1,3),'.m')
% $$$ plot3(xyz(sind,6,1),xyz(sind,6,2),xyz(sind,6,3),'.r')
% $$$ plot3(xyz(sind,7,1),xyz(sind,7,2),xyz(sind,7,3),'.g')
% $$$ plot3(xyz(sind,8,1),xyz(sind,8,2),xyz(sind,8,3),'.k')
% $$$ plot3(xyz(sind,8,1),xyz(sind,8,2),xyz(sind,8,3),'.k')
% $$$ plot3(nmark(sind,1),nmark(sind,2),nmark(sind,3),'*r')
% $$$ daspect([1,1,1])
% $$$ 
% $$$ 
% $$$ 
% $$$ scatter3(xyz(ind,1,1),xyz(ind,1,2),xyz(ind,1,3),30,'m')
% $$$ plot3(xyz(ind,6,1),xyz(ind,6,2),xyz(ind,6,3),'.r')
% $$$ plot3(xyz(ind,4,1),xyz(ind,4,2),xyz(ind,4,3),'.c')
% $$$ plot3(xyz(ind,7,1),xyz(ind,7,2),xyz(ind,7,3),'.g')
% $$$ plot3(xyz(ind,8,1),xyz(ind,8,2),xyz(ind,8,3),'.k')


cxyz = cat(2,xyz(ind,[1],:),bvecs,xyz(ind,[6,7,8],:));
c = [0,1,0;...
     0.5,0,0.5;...
     1,0,0;...
     0,0,1;...
     0,0,1;...
     0,0,1];



cind = 1000;
figure
hax = axes;
marker_options = {};
marker_options.style ='o' ;
marker_options.size = 8 ;
marker_options.erase = 'none'; %{none|xor|background}
markers = repmat({0},1,size(cxyz,2));

for l=1:size(cxyz,2),
    for b=l:size(cxyz,2),
        sticks{l,b} =animatedline([cxyz(cind(end),l,1),...
                                    cxyz(cind(end),b,1)],...
                                [cxyz(cind(end),l,2),...
                                    cxyz(cind(end),b,2)],...
                                [cxyz(cind(end),l,3),...
                                    cxyz(cind(end),b,3)]);
       set(sticks{l,b},'Parent',hax,...
                     'LineWidth', 2, ... stick_options.width,
                     'Visible','on',...
                     'Parent',hax);
        
    end
end

% $$$ c = [0,1,0;...
% $$$      0,0,1;...
% $$$      0,0,1;...
% $$$      1,0,0;...
% $$$      0.5,0,0.5;...
% $$$      0,1,0;...
% $$$      0,0,1;...
% $$$      0,0,1;...
% $$$      1,0,0];

for l=1:size(cxyz,2),
    markers{l} =animatedline([cxyz(cind(end),l,1),cxyz(cind(end),l,1)],...
                             [cxyz(cind(end),l,2),cxyz(cind(end),l,2)],...
                             [cxyz(cind(end),l,3),cxyz(cind(end),l,3)]);
    set(markers{l},'Parent',hax, ...
                   'Marker'         , marker_options.style, ...
                   'MarkerEdgeColor', c(l,:), ...
                   'MarkerSize'     , marker_options.size, ...
                   'MarkerFaceColor', c(l,:), ...
                   'Visible','on');
    
end    


xlim([-0.6, 0.6])
ylim([-0.5, 0.5])
zlim([0, 0.6])
daspect([1,1,1])

for cind = 1:10:10000
    for l=1:size(cxyz,2),
        for b=l:size(cxyz,2),

            clearpoints(sticks{l,b});
            addpoints(sticks{l,b},[cxyz(cind(end),l,1),...
                                   cxyz(cind(end),b,1)],...
                                  [cxyz(cind(end),l,2),...
                                   cxyz(cind(end),b,2)],...
                                  [cxyz(cind(end),l,3),...
                                   cxyz(cind(end),b,3)]);
            
        end
        clearpoints(markers{l});
        addpoints(markers{l},[cxyz(cind(end),l,1),cxyz(cind(end),l,1)],...
                             [cxyz(cind(end),l,2),cxyz(cind(end),l,2)],...
                             [cxyz(cind(end),l,3),cxyz(cind(end),l,3)]);
    
    end
    drawnow
    pause(0.001);
end