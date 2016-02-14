function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_tsne_rev15(Trial,varargin)
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


fet.data = [fxyz(:,{'pelvis_root','spine_middle','spine_upper','hcom'},3),...  1,2,3,4
            fvelxy(:,{'spine_lower','spine_middle'}),...                   5,6
            fang(:,'spine_lower','pelvis_root',2),...                      7            
            fang(:,'pelvis_root','spine_middle',2),...                     8            
            fang(:,'spine_middle','spine_upper',2),...                     9            
            zv.data,...                                                    10
            sv.data,...                                                    11
            av.data,...                                                    12
            zav.data,...                                                   13
            fvelz(:,{'spine_upper'}),...                                   14
            abs(circ_dist(circ_mean([fang(:,1,4,1),fang(:,1,3,1),fang(:,1,2,1)],[],2),...15
                          fang(:,4,'hcom',1))),...
            fang(:,1,4,3).*cos(fang(:,1,4,2))./prctile(fang(:,1,4,3).*cos(fang(:,1,4,2)),98)...16
            ];
fet.data(isinf(fet(:))) = 0;


if normalize,
    fet.unity;
end


%%% REDO LABELS AND TITLES
% This should be converted to a model elements
featureTitles = {};
featureDesc = {};
if nargout>1,
    %% Feature tags and definitions
    %lower spine speed
    % 1.
    %  '$\displaystyle\frac{dz{BL}}{dt}$'
    featureTitles(end+1) = {'Z BP'};
    %featureTitles(end+1) = {'Z_{BP}'};
    featureDesc(end+1) = {'1 Hz low pass filtered height of the pelvic maker'};

    % 2 
    featureTitles(end+1) = {'Z BM'};
    %featureTitles(end+1) = {'Z_{BM}'};
    featureDesc(end+1) = {'1 Hz low pass filtered height of the middle spine maker'};

    % 3
    featureTitles(end+1) = {'Z BU'};
    %featureTitles(end+1) = {'Z_{BU}'};
    featureDesc(end+1) = {'1 Hz low pass filtered height of the middle spine maker'};

    % 4
    featureTitles(end+1) = {'Z HC'};
    %featureTitles(end+1) = {'Z_{HC}'};
    featureDesc(end+1) = {['1 Hz low pass filtered height of the ' ...
                        'head''s center of mass']};

    % 5. 
    featureTitles(end+1) = {'XY Speed BL'};
    %featureTitles(end+1) = {'d(XY_{BL})/dt'};
    %featureTitles(end+1) = {'$log10(\displaystyle\frac{dZ_{BL}}{dt}) \quad log10(mm/s)$'};
    featureDesc(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                    'the spine lower maker']};
    % 6.
    featureTitles(end+1) = {'XY Speed BM'};
    %featureTitles(end+1) = {'XY Speed_{BM}'};
    %featureTitles(end+1) = {'$log10(\displaystyle\frac{dZ_{BM}}{dt}) \quad log10(mm/s)$'};    
    featureDesc(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                    'the middle spine marker']};
    % 7.
    featureTitles(end+1) = {'PPC traj yaw'};    
    %featureTitles(end+1) = {'PPC_{traj yaw}'};
    featureDesc(end+1) = {['1 Hz lowpass filtered Pair-wise Phase Consisistency(PPC) ',...
                           'of the yaw of trajectories of all makers along the ',...
                           'rostro-caudal axis']};
    % 8.
    featureTitles(end+1) = {'Pitch BMBU'};    
    %featureTitles(end+1) = {'Pitch_{BMBU}'};
    %featureTitles(end+1) = {'$\psi_{BMBU} \quad rad$'};
    featureDesc(end+1) = {['Pitch of body middle to body upper relative to xy ' ...
                           'plane']};
    % 9.
    featureTitles(end+1) = {'Pitch BLBM'};    
    %featureTitles(end+1) = {'Pitch_{BLBM}'};
    %featureTitles(end+1) = {'$\psi_{BLBM} \quad rad$'};
    featureDesc(end+1) = {['Pitch of body lower to body middle relative to xy ' ...
                           'plane']};
    % 10.
    featureTitles(end+1) = {'Pitch HBHF'};    
    %featureTitles(end+1) = {'Pitch_{HBHF}'};
    %featureTitles(end+1) = {'$\psi_{HBHF} \quad rad$'};
    featureDesc(end+1) = {['Pitch of head_back to head_front relative to xy ' ...
                           'plane']};

    % 11.
    featureTitles(end+1) = {'|BL p-> BLBU|'};
    featureDesc(end+1) = {['Magnitude of the projection of lower spine trajectory  ' ...
                           'onto the vecor of lower spine to upper spine']};
    % 12.
    featureTitles(end+1) = {'Spine Sinuosity'};
    featureDesc(end+1) = {['Length of the spine divided by the ' ...
                        'distance beween the endpoints']};

    % 13.
    featureTitles(end+1) = {'XYZ Dist BLBM'};
    %featureTitles(end+1) = {'XYZ Dist_{BLBM}'};
    %featureTitles(end+1) = {'$r_{BLBM} \quad mm$'};
    featureDesc(end+1) = {['Distance between body lower and body middle.']};


    % 14.
    featureTitles(end+1) = {'CV^-1 \theta spine'};
    featureDesc(end+1) = {['Inverse coefficient of variation of ' ...
                        'spine angular speed']};

    % 11.
    featureTitles(end+1) = {'XYZ Dist BPBU'};
    %featureTitles(end+1) = {'XYZ Dist_{BPBU}'};
    %featureTitles(end+1) = {'$r_{BPBU} \quad mm$'};
    featureDesc(end+1) = {['Distance between body pelvis and body upper.']};

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
    featureTitles(end+1) = {'XYZ Dist BMBU'};
    %featureTitles(end+1) = {'XYZ Dist_{BMBU}'};
    %featureTitles(end+1) = {'$r_{BMBU} \quad mm$'};
    featureDesc(end+1) = {['Distance between body midle and body upper.']};

    % 15.
    featureTitles(end+1) = {'Z Speed BU'};
    featureDesc(end+1) = {['1 Hz low pass filtered speed in the z axis of the ' ...
                    'upper spine marker']};

    % 16.  
    featureTitles(end+1) = {'d(yaw BLBU)/dt'};
    %featureTitles(end+1) = {'d(yaw_{BLBU})/dt'};
    %featureTitles(end+1) = {'$log10(\displaystyle\frac{d\theta_{BLBU}}{dt}) \quad log10(rad/s)$'};
    featureDesc(end+1) = {'Pitch speed of the vector from spine_lower to spine_upper'};

    % 17.
    featureTitles(end+1) = {'XY Dist BLBU'};
    %featureTitles(end+1) = {'XY Dist_{BLBU}'};
    %featureTitles(end+1) = {'$normalized(r_{BPBM}) \quad A.U.$'};        
    featureDesc(end+1) = {['Normalized distance in the xy plane between the upper and ',...
                           'lower spine.']};
    % 18.
    featureTitles(end+1) = {'XYZ Dist BLHC'};
    %featureTitles(end+1) = {'XYZ Dist_{BLHC}'};
    %featureTitles(end+1) = {'$normalized(r_{BPBM}) \quad A.U.$'};    
    featureDesc(end+1) = {['Normalized distance between head and tail.']};

end
