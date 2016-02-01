function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_tsne_rev13(Trial,varargin)
[newSampleRate,normalize] = DefaultArgs(varargin,{15,false},1);

if ischar(Trial),
    Trial = MTATrial(Trial);
elseif iscell(Trial),
    Trial = MTATrial(Trial{:});
end



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
mxyz = Trial.load('xyz');
mxyz.filter('ButFilter',3,20);
tsh = round(.15*mxyz.sampleRate);
afet = mxyz.copy;
afet.data = circshift(mxyz.data,-tsh)-circshift(mxyz.data,tsh);
afet.data = reshape(afet.data,[],3);
afet.data = permute(bsxfun(@dot,permute(reshape(repmat(mxyz(:,4,:)-mxyz(:,1,:),[1,mxyz.size(2),1]),[],3),[2,1]),permute(afet.data,[2,1])),[2,1]);
afet.data = reshape(afet.data,[],mxyz.size(2));
afet.resample(newSampleRate);
zv = afet.copy;
zv.data = log10(abs(afet(:,1))).*sign(afet(:,1));





%% XY speed
fvelxy = xyz.copy;
fvelxy.filter('ButFilter',3,2.4,'low');
fvelxy = fvelxy.vel([],[1,2]);
fvelxy.resample(newSampleRate);
fvelxy.data(fvelxy.data<=10e-5) = 10e-5;
fvelxy.data = log10(fvelxy.data);



%% Z speed
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,1,'low');
fvelz = fxyz.vel([],[3]);
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


%% ANG Vel
sh = round(.2*fang.sampleRate);
ang_b_vel = [abs(circ_dist(circshift(fang(:,1,4,1),-sh),circshift(fang(:,1,4,1),sh)))];
ang_b_vel(ang_b_vel<.1e-6) = 1e-6;
ang_b_vel = log10(ang_b_vel);


%% SS 
sd = sqrt(sum((ss.data-circshift(ss.data,-1,2)).^2,3));
sv = Trial.xyz.copy;
sv.data = sum(sd(:,2:end-1),2)./sd(:,end);
sv.resample(fxyz);


%% AV
sang = [circ_dist(fang(:,1,2,1),fang(:,2,3,1)),...
        circ_dist(fang(:,2,3,1),fang(:,3,4,1)),...
        circ_dist(fang(:,3,4,1),fang(:,4,'hcom',1))];
av = fang.copy;
av.data = abs(sum(sang,2)-circ_mean(sum(sang,2)));


%% ZAV
sh = 1;
sang = [circ_dist(circshift(fang(:,1,2,1),-sh),circshift(fang(:,1,2,1),sh)),...
        circ_dist(circshift(fang(:,1,3,1),-sh),circshift(fang(:,1,3,1),sh)),...
        circ_dist(circshift(fang(:,1,4,1),-sh),circshift(fang(:,1,4,1),sh)),...
        circ_dist(circshift(fang(:,1,5,1),-sh),circshift(fang(:,1,5,1),sh)),...        
        circ_dist(circshift(fang(:,1,7,1),-sh),circshift(fang(:,1,7,1),sh))...                
       ];
zav = Trial.xyz.copy;
zav.data =  log10(mean(abs(sang),2)./(var(sang,[],2)+1));

%% Main feature construction
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'tSNE_Features',mfilename,'t');                  


fet.data = [fxyz(:,{'pelvis_root','spine_middle','spine_upper','hcom'},3),...            
            fvelxy(:,{'spine_lower','spine_middle'}),...
            man.data+.2,...
            fang(:,'spine_middle','spine_upper',2),...                  %Pitch            
            fang(:,'spine_lower','spine_upper',2),...                  %Pitch            
            fang(:,'head_back','head_front',2),...
            zv.data,...
            sv.data,...
            av.data,...
            zav.data,...
            fvelz(:,{'spine_upper'}),...
            ang_b_vel,...
            abs(circ_dist(circ_mean([fang(:,1,4,1),fang(:,1,3,1),fang(:,1,2,1)],[],2),...
                          fang(:,4,'hcom',1))),...
            fang(:,1,4,3).*cos(fang(:,1,4,2))./prctile(fang(:,1,4,3).*cos(fang(:,1,4,2)),98)...
            ];%rhm.data
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

    featureTitles(end+1) = {'Height_{BM}'};
    featureDesc(end+1) = {'1 Hz low pass filtered height of the upper spine maker'};
    % 2.
    featureTitles(end+1) = {'XY Dist_{BLBU}'};
    featureDesc(end+1) = {['distance in the xy plane between the upper and ',...
                           'lower spine.']};
    % 6.
    featureTitles(end+1) = {'XY Speed_{BL}'};
    featureDesc(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                    'the spine lower maker']};
    % 7.
    featureTitles(end+1) = {'XY Speed_{BU}'};
    featureDesc(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                    'the spine upper maker']};
    % 8.
    featureTitles(end+1) = {'XY Speed_{HC}'};
    featureDesc(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                    'the center of mass of the head']};
    % 9.
    featureTitles(end+1) = {'Z Speed(flp1Hz) of Upper Spine'};
    featureDesc(end+1) = {['1 Hz low pass filtered speed in the z axis of the ' ...
                    'upper spine marker']};
    % 10.
    featureTitles(end+1) = {'PPC_{traj yaw}'};
    featureDesc(end+1) = {['1 Hz lowpass filtered Pair-wise Phase Consisistency(PPC) ',...
                           'of the yaw of trajectories of all makers along the ',...
                           'rostro-caudal axis']};
    % 11.
    featureTitles(end+1) = {'bfet'};
    featureDesc(end+1) = {['Magnitude of the projection of lower spine trajectory  ' ...
                    'onto the vecor of lower spine to upper spine']};
    % 12.
    featureTitles(end+1) = {'zfet'};
    featureDesc(end+1) = {['Spectral power of low frequency(<7Hz) oscillations along the ' ...
                    'z-axis of the lower spine marker']};
    % 13.
    featureTitles(end+1) = {'Pitch_{BMBU}'};
    featureDesc(end+1) = {['Pitch of spine_middle to spine_upper relative to xy ' ...
                    'plane']};
    % 14.
    featureTitles(end+1) = {'Pitch_{BUHB}'};
    featureDesc(end+1) = {['Pitch of spine_upper to head_back relative to xy ' ...
                    'plane']};
    % 15.
    featureTitles(end+1) = {'Pitch_{HBHF}'};
    featureDesc(end+1) = {['Pitch of head_back to head_front relative to xy ' ...
                    'plane']};
    % 16.
    featureTitles(end+1) = {'d(pitch_{BMBU})/dt'};
    featureDesc(end+1) = {'Pitch speed of the vector from spine_middle to spine_upper'};
    % 17.
    featureTitles(end+1) = {'d(yaw_{BLBU})/dt'};
    featureDesc(end+1) = {'Pitch speed of the vector from spine_lower to spine_upper'};
    % 18.
    featureTitles(end+1) = {'d(yaw_{BMHF})/dt'};
    featureDesc(end+1) = {'Yaw speed of the vector from spine_middle to head_front'};
    % 19.
    featureTitles(end+1) = {'|yaw|_{BL,BM,HB}'};
    featureDesc(end+1) = {'angle between head, spine and tail'};
    % 20.
    featureTitles(end+1) = {'XYZ Dist_{BLHC}'};
    featureDesc(end+1) = {['Distance between head and tail.']};

end

