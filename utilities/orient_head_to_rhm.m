function orient_head_to_rhm(filename,varargin)


%filename = '/storage/nickdg/data/VR Engagement/justin_VRObj_ratdata.mat';
%outputFileName = '/storage/gravio/ownCloud/Shared/VR_Methods/matlab/justin_VRObj_ratdata_CorrectedOrientation.mat';
%filename = '/storage/nickdg/data/Wall Avoidance/justin_WallAvoid_ratdata.mat';
%outputFileName = '/storage/gravio/ownCloud/Shared/VR_Methods/matlab/justin_WallAvoid_ratdata_CorrectedOrientation_spwFiltered.mat';
%outputFileName = '/storage/gravio/ownCloud/Shared/VR_Methods/matlab/justin_WallAvoid_ratdata_CorrectedOrientation_spwFilterednsr.mat';
%filename ='/storage/gravio/ownCloud/Shared/VR_Methods/matlab/justin_WallAvoid_ratdata_CorrectedOrientation.mat';
%outPfilename ='/storage/gravio/ownCloud/Shared/VR_Methods/matlab/justin_WallAvoid_ratdata_CorrectedOrientation.mat';
%filename =
%'/storage/gravio/ownCloud/Shared/VR_Methods/matlab/justin_VRObj_ratdata_CorrectedOrientation.mat';


[dir,filename,ext] = fileparts(filename);
load(fullfile(dir,[filename,ext]));


if ~isempty(varargin)
    outputFileName = varargin{1};
    [odir,ofilename,oext] = fileparts(outputFileName);    
else
    odir = dir;
    ofilename = filename;
    oext = ext;
end

dispFigs = true;

if dispFigs,
    set(0,'defaultAxesFontSize',10,...
          'defaultTextFontSize',10)
    mkdir(fullfile(odir,ofilename));
end
%data(:,1): time
% rigidbody center of mass
%data(:,2): Y    Nick's X
%data(:,3): Z    Nick's Y
%data(:,4): X    Nick's Z
% euler angles of rigid body relative to room cooridate system
%data(:,5): Roll
%data(:,6): Pitch
%data(:,7): Yaw
% Some orientation vector nick made
%data(:,8): Y
%data(:,9): Z
%data(:,10): X
% Session Id
%data(:,11): session id


% euler ang vector
eang = deg2rad(data(:,[5,6,7]));

yaw = eang(:,1);
pitch = eang(:,2);
roll = eang(:,3);

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
sessionIds = unique(data(:,11));

%xyzSampleRate = 150;
xyzSampleRate = 240;
% store XYZ rigid body center of mass in meters
% X -> long axis of maze
% Y -> short axis of maze
% Z -> Height
% xyz(time,marker,dimension)
xyz(:,1,:) = data(:,[4,2,3]);
% Add a low passed (0.8 Hz) version of center of mass with slight time shift
xyz(:,2,:) = circshift(ButFilter(xyz(:,1,:),3,[.8]./(xyzSampleRate/2),'low'),round(xyzSampleRate/4));
% Add a low passed (0.2 Hz) version of center of mass with large time shift
xyz(:,3,:) = circshift(ButFilter(xyz(:,1,:),3,[.2]./(0.5*xyzSampleRate),'low'),round(xyzSampleRate));
% Add orientation vector + center of mass reduce in length
xyz(:,4,:) = 0.05.*data(:,[10,8,9])+data(:,[4,2,3]);
% Add a low passed (2 Hz) version of center of mass with no time shift
xyz(:,5,:) = ButFilter(xyz(:,1,:),3,[2]./(xyzSampleRate/2),'low');
% Add markers surrounding the center of mass to create a head basis
xyz(:,6,:) =  sq(evec(:,1,:))+data(:,[4,2,3]);
xyz(:,7,:) =  sq(evec(:,2,:))+data(:,[4,2,3]);
xyz(:,8,:) =  sq(evec(:,3,:))+data(:,[4,2,3]);


% create low passed (2.4 Hz) speed in xy plane
markers = [1,5];
vxy = sqrt(sum([zeros([1,numel(markers),2]);diff(xyz(:,markers,[1,2]))].^2,3)).*xyzSampleRate.*100;
vxy(vxy(:,2)<1e-3,:) = 1e-3;
vxy = log10(vxy);

% Find Best orientation (good luck)

for ses = 1:numel(sessionIds),

    ind = data(:,11)==sessionIds(ses);

    %% Display the Original basis for the head
    if dispFigs,
        hfig = figure(20392039);
        set(hfig,'units','centimeters')
        set(hfig,'Position',[11 0 32 25])
        set(hfig,'PaperPositionMode','auto');
        clf;
        subplot(2,2,1);
        hold('on');
        moffset = 1200;
        sind = find(ind,1,'first')+moffset;
        plot3(xyz(sind,1,1),xyz(sind,1,2),xyz(sind,1,3),'.m')
        plot3(xyz(sind,6,1),xyz(sind,6,2),xyz(sind,6,3),'.r')
        plot3(xyz(sind,7,1),xyz(sind,7,2),xyz(sind,7,3),'.g')
        plot3(xyz(sind,8,1),xyz(sind,8,2),xyz(sind,8,3),'.k')
        title('Vector Search Space')
    end


    
    %% Rotate a basis around the unit sphere
    spw = [];
    tots = 0:5:360;
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
            ovec = nvec(:,2,:);
            ovec = cross(nvec(:,2,:),nvec(:,4,:));
            ovec = bsxfun(@rdivide,ovec,sqrt(sum(ovec.^2,3))).*blen;

            % Create new marker
            tmark = permute(sum(head_rotMat.*permute(reshape(ovec(:,j(ones(3,1),:)),[size(head_norm,1),3,3]),[2,3,1]),2),[3,1,2])+data(ind,[4,2,3]);
            sxyz = cat(2,xyz(ind,5,:),permute(tmark,[1,3,2]));
            
            %% Calculate RHM power for vector
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

             if dispFigs,
                 plot3(tmark(moffset,1),tmark(moffset,2),tmark(moffset,3),'*r');
             end
        end
    end


    
    
    %% Automatic orientation detection
    % spw is the map where the peaks represents the front and back of
    % the head. The indicies (d,t) are the two angle thingies you need 
    % to get the final orientation vector. 

    nspw = spw(nniz(spw),:);
    fspw = imgaussfilt(repmat(nspw,3,3),2);
    fspw = fspw(size(nspw,1)+1:size(nspw,1)*2,size(nspw,2)+1:size(nspw,2)*2);
    
    if dispFigs,
        subplot(222);
        imagesc(nspw');
        title('Rhythmic(6-13Hz) mean PSD')
        xlabel('\phi (rad)')
        ylabel('\theta (rad)');
    end
    [mins,minv] = LocalMinimaN(-fspw,-mean(fspw(:)),10);

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

        tmark = permute(permute(sum(head_rotMat.*permute(reshape(ovec(:,j(ones(3,1),:)),[size(head_norm,1),3,3]),[2,3,1]),2),[3,1,2])+data(ind,[4,2,3]),[1,3,2]);

        bvecs = cat(2,bvecs,tmark);
    end


    %% Calculate RHM power    
    cxyz = cat(2,xyz(ind,[1,2,3],:),bvecs,xyz(ind,[6,7,8],:));

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
        if dispFigs&&m<=2,
            subplot(2,2,2+m);
            vind = nniz(vxy(ind,2));
            tvxy = vxy(ind,2);
            hist2([circ_dist(ang(vind,2,3,1),ang(vind,1,m+3,1)),tvxy(vind)],40,linspace(-.5,2,40));
            title(['Vector ',num2str(m)]);
            ylabel('log10 speed (cm/s)');
            xlabel('Angle between Movement vs Vector direction (rad)');
            %rose(circ_dist(ang(vind,2,1,1),ang(vind,1,m+3,1)),100);
        end
    end

    m = min(abs(mang))==abs(mang);    
    if round(abs(mang(m)),0)==3,
        revOri = -1;
    else
        revOri = 1;
    end

    dout = permute(bvecs(:,min(abs(mang))==abs(mang),:),[1,3,2])-data(ind,[4,2,3]);
    data(ind,8:10) = revOri.*dout(:,[2,3,1]);

    if dispFigs,
        suptitle([filename,' SessionId: ',num2str(ses)])
        print(hfig,'-depsc2',fullfile(odir,ofilename,[ofilename,'-sessionId-',num2str(ses),'.eps']));
        print(hfig,'-dpng',fullfile(odir,ofilename,[ofilename,'-sessionId-',num2str(ses),'.png']));
        delete(hfig)
    end

    disp(['Completed session: ' num2str(ses)])
end




%% Save the file
if ~isempty(outputFileName),
    save(outputFileName,'data');
else
    save(fullfile(dir,[filename,'_CorrectedOrientation',ext]),'data');
end



%% Double Check which direction the head is pointed relative to
%  the direction of movement
if dispFigs,
    hfig = figure(230302030);
    set(hfig,'units','centimeters')
    set(hfig,'Position',[11 0 32 25])
    set(hfig,'PaperPositionMode','auto');

    sfac = factor(numel(sessionIds));

    while numel(sfac)>2,
        sfac = sort(sfac,'descend');
        sfac(end-1) = sfac(end)*sfac(end-1);
        sfac(end) = [];
    end
    if numel(sfac)==1
        sfac = [sfac,1];
    end

    
    
    for ses = 1:numel(sessionIds),

        ind = data(:,11)==sessionIds(ses);
        cxyz = cat(2,xyz(ind,[1,2,3],:),permute(data(ind,[10,8,9])+data(ind,[4,2,3]),[1,3,2]));
        %cxyz = cat(2,xyz(ind,[1,2],:),bvecs(:,1,:));

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


        vind = nniz(vxy(ind,2));
        tvxy = vxy(ind,2);

        subplot(sfac(2),sfac(1),ses)
        hist2([circ_dist(ang(vind,3,2,1),ang(vind,1,4,1)),tvxy(vind)],40,linspace(-.5,2,40));
        title(['session: ',num2str(ses)])
    end
    print(hfig,'-depsc2',fullfile(odir,ofilename,[ofilename,'summary.eps']));
    print(hfig,'-dpng',fullfile(odir,ofilename,[ofilename,'summary.png']));
end 





% $$$ cxyz = cat(2,xyz(ind,[1],:),bvecs,xyz(ind,[6,7,8],:));
% $$$ c = [0,1,0;...
% $$$      0.5,0,0.5;...
% $$$      1,0,0;...
% $$$      0,0,1;...
% $$$      0,0,1;...
% $$$      0,0,1];
% $$$ 
% $$$ 
% $$$ 
% $$$ cind = 1000;
% $$$ figure
% $$$ hax = axes;
% $$$ marker_options = {};
% $$$ marker_options.style ='o' ;
% $$$ marker_options.size = 8 ;
% $$$ marker_options.erase = 'none'; %{none|xor|background}
% $$$ markers = repmat({0},1,size(cxyz,2));
% $$$ 
% $$$ for l=1:size(cxyz,2),
% $$$     for b=l:size(cxyz,2),
% $$$         sticks{l,b} =animatedline([cxyz(cind(end),l,1),...
% $$$                                     cxyz(cind(end),b,1)],...
% $$$                                 [cxyz(cind(end),l,2),...
% $$$                                     cxyz(cind(end),b,2)],...
% $$$                                 [cxyz(cind(end),l,3),...
% $$$                                     cxyz(cind(end),b,3)]);
% $$$        set(sticks{l,b},'Parent',hax,...
% $$$                      'LineWidth', 2, ... stick_options.width,
% $$$                      'Visible','on',...
% $$$                      'Parent',hax);
% $$$         
% $$$     end
% $$$ end
% $$$ 
% $$$ % $$$ c = [0,1,0;...
% $$$ % $$$      0,0,1;...
% $$$ % $$$      0,0,1;...
% $$$ % $$$      1,0,0;...
% $$$ % $$$      0.5,0,0.5;...
% $$$ % $$$      0,1,0;...
% $$$ % $$$      0,0,1;...
% $$$ % $$$      0,0,1;...
% $$$ % $$$      1,0,0];
% $$$ 
% $$$ for l=1:size(cxyz,2),
% $$$     markers{l} =animatedline([cxyz(cind(end),l,1),cxyz(cind(end),l,1)],...
% $$$                              [cxyz(cind(end),l,2),cxyz(cind(end),l,2)],...
% $$$                              [cxyz(cind(end),l,3),cxyz(cind(end),l,3)]);
% $$$     set(markers{l},'Parent',hax, ...
% $$$                    'Marker'         , marker_options.style, ...
% $$$                    'MarkerEdgeColor', c(l,:), ...
% $$$                    'MarkerSize'     , marker_options.size, ...
% $$$                    'MarkerFaceColor', c(l,:), ...
% $$$                    'Visible','on');
% $$$     
% $$$ end    
% $$$ 
% $$$ 
% $$$ xlim([-0.6, 0.6])
% $$$ ylim([-0.5, 0.5])
% $$$ zlim([0, 0.6])
% $$$ daspect([1,1,1])
% $$$ 
% $$$ for cind = 1:10:10000
% $$$     for l=1:size(cxyz,2),
% $$$         for b=l:size(cxyz,2),
% $$$ 
% $$$             clearpoints(sticks{l,b});
% $$$             addpoints(sticks{l,b},[cxyz(cind(end),l,1),...
% $$$                                    cxyz(cind(end),b,1)],...
% $$$                                   [cxyz(cind(end),l,2),...
% $$$                                    cxyz(cind(end),b,2)],...
% $$$                                   [cxyz(cind(end),l,3),...
% $$$                                    cxyz(cind(end),b,3)]);
% $$$             
% $$$         end
% $$$         clearpoints(markers{l});
% $$$         addpoints(markers{l},[cxyz(cind(end),l,1),cxyz(cind(end),l,1)],...
% $$$                              [cxyz(cind(end),l,2),cxyz(cind(end),l,2)],...
% $$$                              [cxyz(cind(end),l,3),cxyz(cind(end),l,3)]);
% $$$     
% $$$     end
% $$$     drawnow
% $$$     pause(0.001);
% $$$ end
