function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_tsne_rev3(Trial,varargin)
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

% $$$ fvelxy = xyz.vel([],[1,2]);
% $$$ fvelxy.resample(newSampleRate);
% $$$ fvelxy.filter('ButFilter',3,2.4,'low');

fvelxy = xyz.copy;
fvelxy.filter('ButFilter',3,2.4,'low');
fvelxy = fvelxy.vel([],[1,2]);
fvelxy.resample(newSampleRate);
fvelxy.data(fvelxy.data<=10e-5) = 10e-5;
fvelxy.data = log10(fvelxy.data);
%fvelxy.data(fvelxy.data<0) = 0;
%fvelxy.data = log10(fvelxy.data);


%% Z speed
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,1,'low');

fvelz = fxyz.vel([],[3]);
fvelz.resample(newSampleRate);
fvelz.filter('ButFilter',3,2.4,'low');
fvelz.data(fvelz.data<=10e-5) = 10e-5;
fvelz.data = log10(fvelz.data);
%fvelz.data(fvelz.data<0) = 0;
%fvelz.data = log10(fvelz.data+1);


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



%% ANG Vel
ang_b_vel = [abs(circ_dist(circshift(fang(:,3,4,2),-1),circshift(fang(:,3,4,2),1))),...
             abs(circ_dist(circshift(fang(:,1,4,1),-1),circshift(fang(:,1,4,1),1))),...
             abs(circ_dist(circshift(fang(:,3,7,1),-1),circshift(fang(:,3,7,1),1)))];
ang_b_vel(ang_b_vel<.1e-6) = 1e-6;
ang_b_vel = log10(ang_b_vel);

% New zfet
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


fet.data = [fxyz(:,{'spine_lower','spine_middle','spine_upper','head_front'},3),...
            fang(:,1,4,3).*cos(fang(:,1,4,2)),...                       %Intermarker Distance
            fvelxy(:,{'spine_lower','spine_upper','head_front'}),...
            fvelz(:,{'spine_upper'}),...
            man.data+.2,...
            bfet.data,...
            zfet.data,...
            fang(:,'spine_middle','spine_upper',2),...                  %Pitch
            fang(:,'spine_lower','head_back',2),...                     %Pitch
            fang(:,'spine_upper','head_front',2),...                      %Pitch
            ang_b_vel,...
            abs(circ_dist(fang(:,1,3,1),circ_mean([fang(:,1,5,1),fang(:,1,7,1)],[],2)))];%rhm.data
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
