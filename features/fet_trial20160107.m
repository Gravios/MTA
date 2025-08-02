function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_trial20160107(Trial,varargin)
[newSampleRate,normalize] = DefaultArgs(varargin,{15,false},1);

if ischar(Trial),
    Trial = MTATrial(Trial);
elseif iscell(Trial),
    Trial = MTATrial(Trial{:});
end


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



xyz.addMarker('bcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'spine_lower','pelvis_root','spine_middle'})));

xyz.addMarker('hcom',[.7,0,.7],{},...
              xyz.com(xyz.model.rb({'head_back','head_left','head_front','head_right'})));


%% XY speed
fvelxy = xyz.copy;
fvelxy.filter('ButFilter',3,2.4,'low');
fvelxy = fvelxy.vel([],[1,2]);
fvelxy.resample(newSampleRate);
fvelxy.data(fvelxy.data<=10e-5) = 10e-5;
fvelxy.data = log10(fvelxy.data);


%% Z speed
fvelz = xyz.copy;
fvelz.filter('ButFilter',3,1,'low');
fvelz = fvelz.vel([],[3]);
fvelz.resample(newSampleRate);
fvelz.filter('ButFilter',3,2.4,'low');
fvelz.data(fvelz.data<=10e-5) = 10e-5;
fvelz.data = log10(fvelz.data);


%% FANG inter marker angles based on filtered xyz
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,1.5,'low');
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


ang_b_vel = [abs(circ_dist(circshift(fang(:,1,4,1),-1),circshift(fang(:,1,4,1),1))),...
             abs(circ_dist(circshift(fang(:,3,'hcom',1),-1),circshift(fang(:,3,'hcom',1),1)))];
tfet = MTADxyz('data',bsxfun(@rdivide,ang_b_vel,10.^fang(:,3,4,2)),'sampleRate',fang.sampleRate);
tfet.filter('ButFilter',3,2,'low');
tfet.data = mean(tfet.data,2);
tfet.data(tfet.data<1e-6) = 1e-6;
tfet.data = log10(tfet.data);


%% Main feature construction
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'tSNE_Features',mfilename,'t');                  


fet.data = [fxyz(:,{'spine_lower','spine_middle','hcom'},3),...
            fvelxy(:,{'spine_lower','spine_middle','hcom'}),...
            [0;diff(fang(:,3,4,2))],...            
            man.data+.2,...
            tfet.data,...
            fang(:,'spine_lower','hcom',3),...   
            fang(:,'spine_lower','spine_middle',2),...
            fang(:,'spine_middle','hcom',2),...
            abs(circ_dist(fang(:,1,3,1),fang(:,1,'hcom',1)))];
fet.data(isinf(fet(:))) = 0;

%fet.data = cat(2,circshift(fet.data,round(fet.sampleRate.*0.25)),fet.data,circshift(fet.data,-round(fet.sampleRate.*0.25)));

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
    %
    featureTitles(end+1) = {'Height_{BL}'};
    featureDesc(end+1) = {'1 Hz low pass filtered height of the lower spine maker'};
    % 2.    
    featureTitles(end+1) = {'Height_{BU}'};
    featureDesc(end+1) = {'1 Hz low pass filtered height of the upper spine maker'};
    % 3.
    featureTitles(end+1) = {'Height_{BU}'};
    featureDesc(end+1) = {'1 Hz low pass filtered height of the upper spine maker'};
    % 4.
    featureTitles(end+1) = {'Height_{HF}'};            
    featureDesc(end+1) = {'1 Hz low pass filtered height of the head front maker'};
    % 5.
    featureTitles(end+1) = {'XY Dist_{BLBU}'};
    featureDesc(end+1) = {['Magnitude of the projection of the vector formed ' ...
                    'by the spine_lower and spine_upper markers']};
    % 6.
    featureTitles(end+1) = {'XY Speed_{BL}'};
    featureDesc(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                    'the spine lower maker']};
    % 7.
    featureTitles(end+1) = {'XY Speed_{BM}'};
    featureDesc(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                    'the spine middle maker']};
    % 7.
    featureTitles(end+1) = {'XY Speed_{BU}'};
    featureDesc(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                    'the spine upper maker']};
    % 8.
    featureTitles(end+1) = {'XY Speed_{HF}'};
    featureDesc(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                    'the head front maker']};
    % 9.
    featureTitles(end+1) = {'Z Speed(flp1Hz) of lower Spine'};
    featureDesc(end+1) = {['1 Hz low pass filtered speed in the z axis of the ' ...
                    'lower spine marker']};
    % 9.
    featureTitles(end+1) = {'Z Speed(flp1Hz) of Middle Spine'};
    featureDesc(end+1) = {['1 Hz low pass filtered speed in the z axis of the ' ...
                    'middle spine marker']};
    % 9.
    featureTitles(end+1) = {'Z Speed(flp1Hz) of Upper Spine'};
    featureDesc(end+1) = {['1 Hz low pass filtered speed in the z axis of the ' ...
                    'upper spine marker']};
    % 9.
    featureTitles(end+1) = {'Z Speed(flp1Hz) of Head Front'};
    featureDesc(end+1) = {['1 Hz low pass filtered speed in the z axis of the ' ...
                    'head front marker']};
    % 10.
    featureTitles(end+1) = {'PPC_{traj yaw}'};
    featureDesc(end+1) = {['1 Hz lowpass filtered Pair-wise Phase Consisistency(PPC) of the yaw of ' ...
                    'trajectories of all makers along the rostro-caudal axis']};
    % 11.
    featureTitles(end+1) = {'bfet'};
    featureDesc(end+1) = {['Magnitude of the projection of lower spine trajectory  ' ...
                    'onto the vecor of lower spine to upper spine']};
    % 12.
    featureTitles(end+1) = {'Pitch_{BMBU}'};
    featureDesc(end+1) = {['Pitch of spine_middle to spine_upper relative to xy ' ...
                    'plane']};
    % 13.
    featureTitles(end+1) = {'Pitch_{BUHB}'};
    featureDesc(end+1) = {['Pitch of spine_upper to head_back relative to xy ' ...
                    'plane']};
    % 14.
    featureTitles(end+1) = {'Pitch_{HBHF}'};
    featureDesc(end+1) = {['Pitch of head_back to head_front relative to xy ' ...
                    'plane']};
    % 15.
    featureTitles(end+1) = {'d(pitch_{BMBU})/dt'};
    featureDesc(end+1) = {'Pitch speed of the vector from spine_middle to spine_upper'};

    % 16.
    featureTitles(end+1) = {'d(yaw_{BLBU})/dt'};
    featureDesc(end+1) = {'Pitch speed of the vector from spine_lower to spine_upper'};
    % 17.
    featureTitles(end+1) = {'d(yaw_{BMHF})/dt'};
    featureDesc(end+1) = {'Yaw speed of the vector from spine_middle to head_front'};
    % 18.
    featureTitles(end+1) = {'|yaw|_{BL,BMH,B}'};
    featureDesc(end+1) = {'angle between head spine and tail'};

end


%% TEsting area %%
% $$$ 
% $$$ figure,
% $$$ plot(abs(sum([circ_dist(fang(:,1,2,1),fang(:,2,3,1)),...
% $$$           circ_dist(fang(:,2,3,1),fang(:,3,4,1)),...
% $$$           circ_dist(fang(:,3,4,1),fang(:,4,7,1))],2)));
% $$$ 
% $$$ 
% $$$ mag = nan([fang.size(1),1]);
% $$$ for a = 1:fang.size(1),
% $$$     mag(a) = PPC([fang(a,1,2,1),fang(a,2,3,1),fang(a,3,4,1),fang(a,4,7,1)]);
% $$$ end
% $$$ 
% $$$ man = fxyz.copy;
% $$$ %man.data = mag;
% $$$ man.data = abs(circ_dist(fang(:,1,3,1),fang(:,1,5,1)));
% $$$ 
% $$$ 
% $$$ 
% $$$ Trial.load('stc','hand_labeled_rev2');
% $$$ 
% $$$ figure,plot(man.data)
% $$$ Lines(Trial.stc{'w',fang.sampleRate}(:),[],'b');
% $$$ Lines(Trial.stc{'n',fang.sampleRate}(:),[],'g');
% $$$ Lines(Trial.stc{'m',fang.sampleRate}(:),[],'m');
% $$$ 
% $$$ eds = linspace(-.2,1,100);
% $$$ eds = linspace(0,3,100);
% $$$ 
% $$$ s = 'm';
% $$$ figure,hold on
% $$$ ind = Trial.stc{['a-' s]};
% $$$ ha = bar(eds,histc(man(ind),eds),'histc');
% $$$ ha.FaceColor = 'c';
% $$$ ha.FaceAlpha = .5;
% $$$ ind = Trial.stc{s};
% $$$ hs = bar(eds,histc(man(ind),eds),'histc');
% $$$ hs.FaceColor = 'r';
% $$$ hs.FaceAlpha = .5;
% $$$ 
% $$$ 
% $$$ 
% $$$ 
