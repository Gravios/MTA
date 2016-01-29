function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_tsne_rev10(Trial,varargin)
[newSampleRate,normalize] = DefaultArgs(varargin,{15,false},1);

if ischar(Trial),
    Trial = MTATrial(Trial);
elseif iscell(Trial),
    Trial = MTATrial(Trial{:});
end

% In case if an old version of MTA constructed the Trial
if isempty(Trial.fet),
    Trial.fet = MTADfet(Trial.spath,...
                        [],...
                        [],...
                        [],...
                        Trial.sync.copy,...
                        Trial.sync.data(1),...
                        []);                  
end


xyz = Trial.load('xyz');

try 
    ss = Trial.load('fet','3dss');
    ss.resample(xyz);
catch
    txyz = xyz.data;
    pnts = zeros([xyz.size(1),105,3]);
    for ind = 1:xyz.size(1),
        try
            pnts(ind,:,:) = fnplt(cscvn(sq(txyz(ind,1:4,:))'))';
        end
    end
    name = '3d spline interpolated spine'; label = '3dss'; key = 's';
    ss = MTADfet.encapsulate(Trial,...
                             pnts,...
                             xyz.sampleRate,...
                             name,label,key);
    ss.updateFilename(Trial);
    ss.save;
end

xyz.data(:,1:4,:) = ss(:,[5,35,65,95],:);

xyz.addMarker('bcom',[.7,0,.7],{},...
    xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})));
xyz.addMarker('hcom',[.7,0,.7],{},...
    xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));





%% TRAJ feature
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,1.5,'low');

bfet = Trial.xyz.copy;
tsh = 20;
afet = Trial.xyz.copy;
afet.data = circshift(fxyz.data,-tsh)-circshift(fxyz.data,tsh);
afet.data = reshape(afet.data,[],3);
afet.data = permute(bsxfun(@dot,permute(reshape(repmat(fxyz(:,4,:)-fxyz(:,1,:),[1,fxyz.size(2),1]),[],3),[2,1]),permute(afet.data,[2,1])),[2,1]);
afet.data = reshape(afet.data,[],fxyz.size(2));
bfet.data = mean(afet(:,1:4),2)./1000.*log10(var(afet(:,1:4),[],2));
bfet.filter('ButFilter',3,2.4,'low');
bfet.resample(newSampleRate);

%% XY speed
fvelxy = xyz.copy;
fvelxy.filter('ButFilter',3,2.4,'low');
fvelxy = fvelxy.vel([],[1,2]);
fvelxy.resample(bfet);
fvelxy.data(fvelxy.data<=10e-5) = 10e-5;
fvelxy.data = log10(fvelxy.data);


%% Z speed
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,1,'low');
fvelz = fxyz.vel([],[3]);
fvelz.resample(bfet);
fvelz.filter('ButFilter',3,2.4,'low');
fvelz.data(fvelz.data<=10e-5) = 10e-5;
fvelz.data = log10(fvelz.data);



%% PPC feature
try
    man = Trial.load('fet','lsppc');
catch err
    gen_fet_lsppc(Trial);    
    man = Trial.load('fet','lsppc');
end
man.filter('ButFilter',3,2,'low');
man.resample(bfet);

%% FANG inter marker angles based on filtered xyz
fxyz = xyz.copy; 
fxyz.filter('ButFilter',3,1.5,'low');
fxyz.resample(bfet);
fang = create(MTADang,Trial,fxyz);


%% ANG Vel
ang_b_vel = [abs(circ_dist(circshift(fang(:,3,4,2),-1),circshift(fang(:,3,4,2),0))),...
             abs(circ_dist(circshift(fang(:,1,4,1),-1),circshift(fang(:,1,4,1),0))),...
             abs(circ_dist(circshift(fang(:,'bcom','hcom',1),-1),circshift(fang(:,'bcom','hcom',1),0)))];
ang_b_vel(ang_b_vel<.1e-6) = 1e-6;
ang_b_vel = log10(ang_b_vel);

%% New zfet
name = 'lower spine Z speed'; label = 'lszs'; key = 'z';
zv = MTADfet.encapsulate(Trial,...
                         diff(xyz(:,1,3)),...
                         xyz.sampleRate,...
                         name,label,key);

dspec = struct('nFFT',2^7,'Fs',zv.sampleRate,...
               'WinLength',2^6,'nOverlap',2^6*.875,...
               'FreqRange',[1,30]);
[ys,fs,ts] = fet_spec(Trial,zv,'mtchglong',false,[],dspec,false);

name = 'zv and wag pw'; label = 'zvwp'; key = 'w';
zfet = MTADfet.encapsulate(Trial,...
                         log10(nanmedian(ys(:,1:fs<7,1,1),2)),...
                         ys.sampleRate,...
                         name,label,key);
zfet.data(~nniz(zfet),:)=-20;
zfet.resample(bfet);



%% Main feature construction
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'tSNE_Features',mfilename,'t');                  


fet.data = [fxyz(:,{'spine_middle'},3),...                 1
            fvelxy(:,{'spine_lower','bcom','hcom'}),...    2,3,4
            man.data+.2,...                                5
            bfet.data,...                                  6
            fang(:,'spine_middle','spine_upper' ,2),...    7
            fang(:,'spine_lower' ,'spine_middle',2),...    8
            fang(:,'head_back'   ,'head_front'  ,2),...    9 
            fang(:,'spine_lower' ,'spine_middle',3),...    10
            fang(:,'pelvis_root' ,'spine_upper' ,3),...    11
            fang(:,'spine_lower' ,'pelvis_root' ,3),...    12
            fang(:,'pelvis_root' ,'spine_middle',3),...    13
            fang(:,'spine_middle','spine_upper' ,3),...    14      
            ang_b_vel,...                                  15,16,17
            abs(circ_dist(fang(:,1,3,1),fang(:,'bcom','hcom',1))),...  18
            fang(:,1,4,3).*cos(fang(:,1,4,2))./prctile(fang(:,1,4,3).*cos(fang(:,1,4,2)),95),... 19
            fang(:,'spine_lower','hcom',3)./prctile(fang(:,'spine_lower','hcom',3),95)];% 20
fet.data(isinf(fet(:))) = 0;


if normalize,
    fet.unity;
end


% This should be converted to a model elements
featureTitles = {};
featureDesc = {};
if nargout>1,
    %% Feature tags and definitions
    %lower spine speed
    % 1.
    %  '$\displaystyle\frac{dz{BL}}{dt}$'
    featureTitles(end+1) = {'Z BM'};
    %featureTitles(end+1) = {'Z_{BM}'};
    featureDesc(end+1) = {'1 Hz low pass filtered height of the upper spine maker'};

    % 2. 
    featureTitles(end+1) = {'d(XY BL)/dt'};
    %featureTitles(end+1) = {'d(XY_{BL})/dt'};
    %featureTitles(end+1) = {'$log10(\displaystyle\frac{dZ_{BL}}{dt}) \quad log10(mm/s)$'};
    featureDesc(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                    'the spine lower maker']};
    % 3.
    featureTitles(end+1) = {'XY Speed BC'};
    %featureTitles(end+1) = {'XY Speed_{BC}'};
    %featureTitles(end+1) = {'$log10(\displaystyle\frac{dZ_{BC}}{dt}) \quad log10(mm/s)$'};    
    featureDesc(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                    'the body''s center of mass']};
    % 4.
    featureTitles(end+1) = {'XY Speed HC'};
    %featureTitles(end+1) = {'XY Speed_{HC}'};
    %featureTitles(end+1) = {'$log10(\displaystyle\frac{dZ_{HC}}{dt}) \quad log10(mm/s)$'};
    featureDesc(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                    'the head''s center of mass']};
    % 5.
    featureTitles(end+1) = {'PPC traj yaw'};    
    %featureTitles(end+1) = {'PPC_{traj yaw}'};
    featureDesc(end+1) = {['1 Hz lowpass filtered Pair-wise Phase Consisistency(PPC) ',...
                           'of the yaw of trajectories of all makers along the ',...
                           'rostro-caudal axis']};
    % 6.
    featureTitles(end+1) = {'bfet'};
    featureDesc(end+1) = {['Magnitude of the projection of lower spine trajectory  ' ...
                           'onto the vecor of lower spine to upper spine']};
    % 7.
    featureTitles(end+1) = {'Pitch BMBU'};    
    %featureTitles(end+1) = {'Pitch_{BMBU}'};
    %featureTitles(end+1) = {'$\psi_{BMBU} \quad rad$'};
    featureDesc(end+1) = {['Pitch of body middle to body upper relative to xy ' ...
                           'plane']};
    % 8.
    featureTitles(end+1) = {'Pitch BLBM'};    
    %featureTitles(end+1) = {'Pitch_{BLBM}'};
    %featureTitles(end+1) = {'$\psi_{BLBM} \quad rad$'};
    featureDesc(end+1) = {['Pitch of body lower to body middle relative to xy ' ...
                           'plane']};
    % 9.
    featureTitles(end+1) = {'Pitch HBHF'};    
    %featureTitles(end+1) = {'Pitch_{HBHF}'};
    %featureTitles(end+1) = {'$\psi_{HBHF} \quad rad$'};
    featureDesc(end+1) = {['Pitch of head_back to head_front relative to xy ' ...
                           'plane']};
    % 10.
    featureTitles(end+1) = {'XYZ Dist BLBM'};
    %featureTitles(end+1) = {'XYZ Dist_{BLBM}'};
    %featureTitles(end+1) = {'$r_{BLBM} \quad mm$'};
    featureDesc(end+1) = {['Distance between body lower and body middle.']};

    % 11.
    featureTitles(end+1) = {'XYZ Dist BPBM'};
    %featureTitles(end+1) = {'XYZ Dist_{BPBM}'};
    %featureTitles(end+1) = {'$r_{BPBM} \quad mm$'};
    featureDesc(end+1) = {['Distance between body pelvis and body middle.']};

    % 12.
    featureTitles(end+1) = {'XYZ Dist BLBP'};
    %featureTitles(end+1) = {'XYZ Dist_{BLBP}'};
    %featureTitles(end+1) = {'$r_{BLBP} \quad mm$'};
    featureDesc(end+1) = {['Distance between body lower and body pelvis.']};

    % 13.
    featureTitles(end+1) = {'XYZ Dist BPBM'};
    %featureTitles(end+1) = {'XYZ Dist_{BPBM}'};
    %featureTitles(end+1) = {'$r_{BPBM} \quad mm$'};
    featureDesc(end+1) = {['Distance between body pelvis and body middle.']};

    % 14.
    featureTitles(end+1) = {'XYZ Dist BPBM'};
    %featureTitles(end+1) = {'XYZ Dist_{BPBM}'};
    %featureTitles(end+1) = {'$r_{BMBU} \quad mm$'};
    featureDesc(end+1) = {['Distance between body midle and body upper.']};

    % 15.
    featureTitles(end+1) = {'d(pitch BMBU)/dt'};    
    %featureTitles(end+1) = {'d(pitch_{BMBU})/dt'};
    %featureTitles(end+1) = {'$log10(\displaystyle\frac{d\psi_{BMBU}}{dt}) \quad log10(rad/s)$'};
    featureDesc(end+1) = {'Pitch speed of the vector from spine_middle to spine_upper'};

    % 17.  
    featureTitles(end+1) = {'d(yaw BLBU)/dt'};
    %featureTitles(end+1) = {'d(yaw_{BLBU})/dt'};
    %featureTitles(end+1) = {'$log10(\displaystyle\frac{d\theta_{BLBU}}{dt}) \quad log10(rad/s)$'};
    featureDesc(end+1) = {'Pitch speed of the vector from spine_lower to spine_upper'};

    % 18.
    featureTitles(end+1) = {'d(yaw BCHC)/dt'};    
    %featureTitles(end+1) = {'d(yaw_{BCHC})/dt'};
    %featureTitles(end+1) = {'$log10(\displaystyle\frac{d\theta_{BCBH}}{dt}) \quad log10(rad/s)$'};
    featureDesc(end+1) = {'Yaw speed of the vector from body COM to head COM'};

    % 19.
    featureTitles(end+1) = {'XY Dist BLBU'};
    %featureTitles(end+1) = {'XY Dist_{BLBU}'};
    %featureTitles(end+1) = {'$normalized(r_{BPBM}) \quad A.U.$'};        
    featureDesc(end+1) = {['Normalized distance in the xy plane between the upper and ',...
                           'lower spine.']};
    % 20.
    featureTitles(end+1) = {'XYZ Dist BLHC'};
    %featureTitles(end+1) = {'XYZ Dist_{BLHC}'};
    %featureTitles(end+1) = {'$normalized(r_{BPBM}) \quad A.U.$'};    
    featureDesc(end+1) = {['Normalized distance between head and tail.']};

end

