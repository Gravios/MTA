function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_tsne(Trial,varargin)
[newSampleRate,normalized,Nmean,Nstd] = DefaultArgs(varargin,{15,false,[],[]},1);


xyz = Trial.load('xyz');

if isempty(Trial.fet),
    Trial.fet = MTADfet(Trial.spath,...
                        [],...
                        [],...
                        [],...
                        Trial.sync.copy,...
                        Trial.sync.data(1),...
                        []);                  
end



%% XYZ Filtered
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,1,'low');

%% TRAJ feature
ang = create(MTADang,Trial,xyz);
tsh = 1;
afet = Trial.xyz.copy;
bfet = Trial.xyz.copy;
afet.data = circshift(xyz(:,:,[1,2]),-tsh)-circshift(fxyz(:,:,[1,2]),tsh);
afet.data = reshape(afet.data,[],2);
aft = mat2cell(afet.data,size(afet,1),[1,1]);
[afet.data,bfet.data] = cart2pol(aft{:});
afet.data = reshape(afet.data,[],xyz.size(2));
bfet.data = reshape(bfet.data,[],xyz.size(2));
m = MTADxyz('data',circ_dist(afet(:,1),ang(:,1,4,1)),'sampleRate',Trial.xyz.sampleRate);
m.data = circ_dist(circshift(m.data,-5),circshift(m.data,5));
m.data = [diff(m.data);0];
bfet.resample(newSampleRate);

%% XY speed
fvelxy = xyz.vel([],[1,2]);
fvelxy.resample(newSampleRate);
fvelxy.filter('ButFilter',3,2.4,'low');

%% Z speed
fvelz = fxyz.vel([],[3]);
fvelz.resample(newSampleRate);
fvelz.filter('ButFilter',3,2.4,'low');

%% FANG inter marker angles based on filtered xyz
fxyz.resample(newSampleRate);
fang = create(MTADang,Trial,fxyz);

%% PPC feature
try
    man = Trial.load('fet','lsppc');
catch err
    gen_fet_lsppc(Trial);    
    man = Trial.load('fet','lsppc');
end
man.filter('ButFilter',3,2,'low');
man.resample(fxyz);



%% RHM feature
[rhm,fs] = fet_rhm(Trial,[],'mtchglong',true);
rhm.data = median(rhm(:,fs>6&fs<14),2);
rhm.resample(fxyz);


fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'tSNE_Features','fet_tsne','t');                  


fet.data = [fxyz(:,{'spine_lower','spine_upper','head_front'},3),...
            fvelxy(:,{'spine_lower','spine_upper','head_front'}),....
            fvelz(:,'head_back'),...
            man.data,...
            log10(abs(bfet(:,1)+1)),...
            fang(:,'spine_middle','spine_upper',2),...
            fang(:,'spine_upper','head_back',2),...
            fang(:,'head_back','head_front',2),...            
            fang(:,1,4,3).*cos(fang(:,1,4,2)),...
            abs(circ_dist(circshift(fang(:,3,4,2),-1),circshift(fang(:,3,4,2),1))),...
            abs(circ_dist(circshift(fang(:,1,4,1),-1),circshift(fang(:,1,4,1),1))),...
            abs(circ_dist(circshift(fang(:,3,7,1),-1),circshift(fang(:,3,7,1),1))),...
            rhm.data];
fet.data(isinf(fet(:))) = 0;


if normalized,
    if isempty(Nmean)||isempty(Nstd),
        [~,Nmean,Nstd] = nunity(fet(Trial.stc{'a'},:),@nan);
    end
    [fet.data] = nunity(fet.data,@nan,Nmean,Nstd);
end

featureTitles = {};
featureDesc = {};
if nargout>1,
    %% Feature tags and definitions
    %lower spine speed
    % 1.
    featureTitles(end+1) = {'Height_{BL}'};
    featureDesc(end+1) = {'1 Hz low pass filtered height of the lower spine maker'};
    % 2.    
    featureTitles(end+1) = {'Height_{BU}'};
    featureDesc(end+1) = {'1 Hz low pass filtered height of the upper spine maker'};
    % 3.
    featureTitles(end+1) = {'Height_{HF}'};            
    featureDesc(end+1) = {'1 Hz low pass filtered height of the head front maker'};
    % 4.
    featureTitles(end+1) = {'XY Speed_{BL}'};
    featureDesc(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                    'the spine lower maker']};
    % 5.
    featureTitles(end+1) = {'XY Speed_{BU}'};
    featureDesc(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                    'the spine upper maker']};
    % 6.
    featureTitles(end+1) = {'XY Speed_{HF}'};
    featureDesc(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                    'the head front maker']};
    % 7.
    featureTitles(end+1) = {'Vertical Speed(flp1Hz) of Middle Spine'};
    featureDesc(end+1) = {['1 Hz low pass filtered speed in the z axis of the ' ...
                    'head back marker']};
    % 8.
    featureTitles(end+1) = {'PPC_{traj yaw}'};
    featureDesc(end+1) = {['1 Hz lowpass filtered Pair-wise Phase Consisistency(PPC) of the yaw of ' ...
                    'trajectories of all makers along the rostro-caudal axis']};
    % 9.
    featureTitles(end+1) = {'bfet'};
    featureDesc(end+1) = {['Magnitude of the projection of lower spine trajectory  ' ...
                    'onto the vecor of lower spine to upper spine']};
    % 10.
    featureTitles(end+1) = {'Pitch_{BMBU}'};
    featureDesc(end+1) = {['Pitch of spine_middle to spine_upper relative to xy ' ...
                    'plane']};
    % 11.
    featureTitles(end+1) = {'Pitch_{BUHB}'};
    featureDesc(end+1) = {['Pitch of spine_upper to head_back relative to xy ' ...
                    'plane']};
    % 12.
    featureTitles(end+1) = {'Pitch_{HBHF}'};
    featureDesc(end+1) = {['Pitch of head_back to head_front relative to xy ' ...
                    'plane']};
    % 13.
    featureTitles(end+1) = {'XY Dist_{BLBU}'};
    featureDesc(end+1) = {['Magnitude of the projection of the vector formed ' ...
                    'by the spine_lower and spine_upper markers']};
    % 14.
    featureTitles(end+1) = {'d(pitch_{BMBU})/dt'};
    featureDesc(end+1) = {'Pitch speed of the vector from spine_middle to spine_upper'};

    % 15.
    featureTitles(end+1) = {'d(yaw_{BLBU})/dt'};
    featureDesc(end+1) = {'Pitch speed of the vector from spine_lower to spine_upper'};
    % 16.
    featureTitles(end+1) = {'d(yaw_{BMHF})/dt'};
    featureDesc(end+1) = {'Yaw speed of the vector from spine_middle to head_front'};
    % 17.
    featureTitles(end+1) = {'rhm'};
    featureDesc(end+1) = {'Rhythmic head motion'};

end