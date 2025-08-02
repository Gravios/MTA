
%filename = '/storage/gravio/ownCloud/Shared/VR_Methods/matlab/justin_WallAvoid_ratdata_CorrectedOrientation.mat';
%filename = '/storage/gravio/ownCloud/Shared/VR_Methods/matlab/justin_WallAvoid_ratdata_CorrectedOrientation_spwFiltered.mat';

%filename = '/storage/gravio/ownCloud/Shared/VR_Methods/matlab/justin_VRObj_ratdata_CorrectedOrientation_spwFiltered.mat';
%filename = '/storage/gravio/ownCloud/Shared/VR_Methods/matlab/justin_VRObj_ratdata_CorrectedOrientation_spwFiltered_nsr.mat';


[dir,filename,ext] = fileparts(filename);
load(fullfile(dir,[filename,ext]));




set(0,'defaultAxesFontSize',10,...
      'defaultTextFontSize',10)
mkdir(fullfile(dir,filename));


% Get the session ids (30 total)

sessionIds = unique(data(:,11));

out.sampleRate =[];
out.xyz = [];
out.ang = [];
out.vel = [];
out.rhm.feature = [];
out.rhm.time = [];
out.rhm.pds = {};
out.rhm.fs  = [];
out.rhm.rhmPow = [];
out.rhm.sampleRate =[];
out.rhm.sesId = [];
out.lagVel = [];
out.angVel = [];
out.tanTraj = [];

for s = sessionIds'    
    ind = data(:,11)==s;    
    xyzSampleRate = 1/median(diff(data(ind,1)));
    out.sampleRate = cat(1,out.sampleRate, xyzSampleRate);
    % store XYZ rigid body center of mass in meters
    % X -> long axis of maze
    % Y -> short axis of maze
    % Z -> Height
    % xyz(time,marker,dimension)
    xyz = [];
    xyz(:,1,:) = data(ind,[4,2,3]);
    % Add orientation vector + center of mass reduce in length
    xyz(:,2,:) = data(ind,[10,8,9])+data(ind,[4,2,3]);
    % Add a low passed (0.8 Hz) version of center of mass with slight time shift
    xyz(:,3,:) = circshift(ButFilter(xyz(:,1,:),3,[.8]./(xyzSampleRate/2),'low'),round(xyzSampleRate/4));
    % Add a low passed (0.2 Hz) version of center of mass with large time shift
    xyz(:,4,:) = circshift(ButFilter(xyz(:,1,:),3,[.2]./(0.5*xyzSampleRate),'low'),round(xyzSampleRate));
    % Add a low passed (2 Hz) version of center of mass with no time shift
    xyz(:,5,:) = ButFilter(xyz(:,1,:),3,[2]./(xyzSampleRate/2),'low');
    % Add a low passed (2 Hz) version of front of head with no time shift
    xyz(:,6,:) = ButFilter(xyz(:,2,:),3,[2]./(xyzSampleRate/2),'low');
    % Add a low passed (0.2 Hz) version of center of mass with large time shift
    xyz(:,7,:) = ButFilter(xyz(:,1,:),3,[.2]./(0.5*xyzSampleRate),'low');
    % obj_center
    xyz(:,8,:) = repmat(permute([-0.00289314048956,0.00922027967768,0.199957333333],[1,3,2]),[sum(ind),1,1]);
    % obj_Corner,
    xyz(:,9,:) = repmat(permute([0.514613546693,0.292023191914,0.213493333333],[1,3,2]),[sum(ind),1,1]);
    % obj_Side,
    xyz(:,10,:) = repmat(permute([0.249139216226,-0.274321778325,0.213493333333],[1,3,2]),[sum(ind),1,1]);

    out.xyz = cat(1,out.xyz,xyz);

    % create low passed (2.4 Hz) speed in xy plane
    markers = [1:6];
    vxy = sqrt(sum([zeros([1,numel(markers),2]);diff(xyz(:,markers,[1,2]))].^2,3)).*xyzSampleRate.*100;
    vxy(vxy(:)<1e-3) = 1e-3;
    vxy = log10(vxy);

    out.vel = cat(1,out.vel,vxy);

    %% Calculate RHM
    fxyz = ButFilter(xyz(:,[5,2],:),3,[50]./(xyzSampleRate/2),'low');

    nframe = size(fxyz,1); %nframe: number of frames (time)
    nmar   = size(fxyz,2);  %nmar: number markers 
    ndim   = size(fxyz,3); %ndim: number of spatial dimensions (fxyz)
    j =1:nmar;
    diffMat = permute(cat(4,permute(reshape(repmat(fxyz(:,:,1)',nmar,1)-fxyz(:,j(ones(nmar,1),:),1).',[nmar,nmar,nframe]),[3,1,2]),permute(reshape(repmat(fxyz(:,:,2)',nmar,1)-fxyz(:,j(ones(nmar,1),:),2).',[nmar,nmar,nframe]),[3,1,2]),permute(reshape(repmat(fxyz(:,:,3)',nmar,1)-fxyz(:,j(ones(nmar,1),:),3).',[nmar,nmar,nframe]),[3,1,2])),[1,3,2,4]);
    ang = zeros(size(fxyz,1),size(fxyz,2),size(fxyz,2),3);
    for i=1:size(fxyz,2),
        for j=1:size(fxyz,2),
            if i==j,continue,end                    
            tang =cell(1,3);
            [tang{:}] = cart2sph(diffMat(:,i,j,1),diffMat(:,i,j,2),diffMat(:,i,j,3));
            ang(:,i,j,:) = cell2mat(tang);
        end
    end
    ang(ang(:,1,2,2)~=0,1,1,1)=1;
    clear('diffMat','nframe','nmar','ndim')


    bang = ButFilter(ang(:,1,2,3),3,[2,50]./(xyzSampleRate/2),'bandpass');
    rhm = [0;ButFilter(diff(bang),3,[2,50]./(xyzSampleRate/2),'bandpass')];


    %% Caculate RHM power spectrum

    parspec = struct('nFFT',2^9,...
                     'Fs',  xyzSampleRate,...
                     'WinLength',2^7,...
                     'nOverlap',2^7*.875,...
                     'NW',3,...
                     'Detrend',[],...
                     'nTapers',[],...
                     'FreqRange',[1,30]);

    [ys,fs,ts] = spec('mtchglong',rhm,parspec);


    % Modify time stamps and spec; add padding (0's)
    ts = ts+(parspec.WinLength/2)/xyzSampleRate;
    ssr = 1/diff(ts(1:2));
    pad = round([ts(1),mod(size(rhm,1)-round(parspec.WinLength/2),parspec.WinLength)/xyzSampleRate].*ssr)-[1,0];
    szy = size(ys);
    rhmpsd =cat(1,zeros([pad(1),szy(2:end)]),ys,zeros([pad(2),szy(2:end)]));
    rhmpsdSampleRate = ssr;
    ts = cat(2,([1:pad(1)]./ssr),ts',([pad(1)+size(ts,1)]+[1:pad(2)])./ssr)';

    rhmpow = median(rhmpsd(:,fs>6&fs<12),2);

% $$$ figure,
% $$$ imagesc(ts,fs,log10(rhmpsd)')
% $$$ axis xy
% $$$ colormap jet


% $$$ xeds = linspace(-.600,.600,60);
% $$$ yeds = linspace(-.400,.400,40);
% $$$ drat = xyzSampleRate/rhmpsdSampleRate;
% $$$ dind = round(linspace(1,size(fxyz,1),numel(rhmpow)));
% $$$ cfx = clip(fxyz(dind,1,1),-.6,.6);
% $$$ cfy = clip(fxyz(dind,1,2),-.4,.4);
% $$$ [~,xinds] = histc(cfx,xeds);
% $$$ [~,yinds] = histc(cfy,yeds);


    out.rhm.feature = cat(1,out.rhm.feature,rhm);
    out.rhm.time = cat(1,out.rhm.time,ts);
    out.rhm.pds = cat(1,out.rhm.pds,{rhmpsd});
    out.rhm.fs  = cat(1,out.rhm.fs,fs);
    out.rhm.rhmPow = cat(1,out.rhm.rhmPow,rhmpow);
    out.rhm.sampleRate = cat(1,out.rhm.sampleRate,rhmpsdSampleRate);
    out.rhm.sesId = cat(1,out.rhm.sesId,s.*ones([size(rhmpsd,1),1]));

% $$$ mrhmxy = accumarray([xinds,yinds],log10(abs(rhmpow)),[numel(xeds),numel(yeds)],@nanmedian);
% $$$ srhmxy = accumarray([xinds,yinds],log10(abs(rhmpow)),[numel(xeds),numel(yeds)],@nanstd);
% $$$ figure,imagesc(xeds,yeds,mrhmxy')
% $$$ colormap defaultp
% $$$ axis xy
% $$$ caxis([-11,-8])
% $$$ 
% $$$ 
% $$$ figure,imagesc(xeds,yeds,srhmxy')
% $$$ colormap jet
% $$$ axis xy
% $$$ caxis([1,2])




    % Locomotion feature

    nframe = size(xyz,1); %nframe: number of frames (time)
    nmar   = size(xyz,2);  %nmar: number markers 
    ndim   = size(xyz,3); %ndim: number of spatial dimensions (xyz)
    j =1:nmar;
    diffMat = permute(cat(4,permute(reshape(repmat(xyz(:,:,1)',nmar,1)-xyz(:,j(ones(nmar,1),:),1).',[nmar,nmar,nframe]),[3,1,2]),permute(reshape(repmat(xyz(:,:,2)',nmar,1)-xyz(:,j(ones(nmar,1),:),2).',[nmar,nmar,nframe]),[3,1,2]),permute(reshape(repmat(xyz(:,:,3)',nmar,1)-xyz(:,j(ones(nmar,1),:),3).',[nmar,nmar,nframe]),[3,1,2])),[1,3,2,4]);
    ang = zeros(size(xyz,1),size(xyz,2),size(xyz,2),3);
    for i=1:size(xyz,2),
        for j=1:size(xyz,2),
            if i==j,continue,end                    
            tang =cell(1,3);
            [tang{:}] = cart2sph(diffMat(:,i,j,1),diffMat(:,i,j,2),diffMat(:,i,j,3));
            ang(:,i,j,:) = cell2mat(tang);
        end
    end
    ang(ang(:,1,2,2)~=0,1,1,1)=1;
    clear('diffMat','nframe','nmar','ndim')

    cang = circ_dist(ang(:,5,4,1),ang(:,2,4,1));
    dcang = circ_dist(circshift(cang,10),circshift(cang,-10));
    bang = clip(ang(:,5,4,3)./sqrt(sq(sum(GetSegs(circshift(dcang,15),1:size(dcang),30).^2)))'./10,0,1000);
    bang(bang~=0&~isnan(bang)) = ButFilter(bang(bang~=0&~isnan(bang)),3,[2.4]./(xyzSampleRate.*0.5),'low');

    out.lagVel = cat(1,out.lagVel,bang);


    out.ang = cat(1,out.ang,ang);
    out.angVel = cat(1,out.angVel,circ_dist(ang(:,4,5,1),circshift(ang(:,4,3,1),10)));



    m = 6;
    dtan = [xyz(:,m,1)-circshift(xyz(:,m,1),-round(0.1*xyzSampleRate)),...
            xyz(:,m,2)-circshift(xyz(:,m,2),-round(0.1*xyzSampleRate)),...
            xyz(:,m,3)-circshift(xyz(:,m,3),-round(0.1*xyzSampleRate))];
    bpr = -sum(dtan.*bsxfun(@rdivide,data(ind,[10,8,9]),sqrt(sum(data(ind,[10,8,9]).^2,2))),2);

    out.tanTraj = cat(1,out.tanTraj,bpr);
end


figure,plot(out.xyz(:,1,1),out.xyz(:,1,2),'.')
hold on,plot(out.xyz(1,8,1),out.xyz(1,8,2),'.m')
hold on,plot(out.xyz(1,9,1),out.xyz(1,9,2),'.m')
hold on,plot(out.xyz(1,10,1),out.xyz(1,10,2),'.m')

figure,hold on
plot(out.ang(:,1,8,3))
plot(out.ang(:,1,9,3))
plot(out.ang(:,1,10,3))


obd = [out.ang(:,1,8,3),...
       out.ang(:,1,9,3),...
       out.ang(:,1,10,3)];
obd = [out.ang(:,1,9,3)];

mdob = min(obd,[],2);

mdt =LocalMinimaN(mdob,0.1,120);

obvs = GetSegs(out.vel(:,1),mdt(:,1)-240,240);
figure,imagesc(obvs')

obts = GetSegs(out.tanTraj(:,1),mdt(:,1)-240,240);
figure,imagesc(obts')

obds = GetSegs(mdob(:,1),mdt(:,1)-240,240);
figure,imagesc(log10(obds)')

rhmSampleRate = 1./diff(out.rhm.time(1:2));
obrhms = GetSegs(out.rhm.rhmPow(:,1),...
                 round((mdt(:,1)-1)/120.*rhmSampleRate)-round(rhmSampleRate.*2),...
                 round(rhmSampleRate.*2));
figure,imagesc(log10(obrhms)')

obdirs = GetSegs(circ_dist(out.ang(:,1,2,1),out.ang(:,8,1,1)),mdt(:,1)-240,240);
figure,imagesc(obdirs')
colormap hsv



mtt = MTADxyz( 'data',out.tanTraj(:,1),'sampleRate',out.sampleRate(1));
mrp = MTADxyz( 'data',out.rhm.rhmPow(:,1),'sampleRate',rhmSampleRate);

mtt.resample(mrp);


mtts = mtt.segs(round((mdt(:,1)-1)/120.*rhmSampleRate)-round(rhmSampleRate.*win),round(rhmSampleRate.*win));
mrps = mrp.segs(round((mdt(:,1)-1)/120.*rhmSampleRate)-round(rhmSampleRate.*win),round(rhmSampleRate.*win));


obds = GetSegs(mdob(:,1),mdt(:,1)-240,240);
(sum(obds,2),100)



mtts(:

figure,
ind = nniz(mtts(:))&nniz(mrps(:));
ind = ':';
hist2([mtts(ind),log10(mrps(ind))],linspace([-.03,.03,50]),linspace([-11,-7.5,50]))



figure



obiPer = [mdt(:,1)-240,mdt(:,1)];

