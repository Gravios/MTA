function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_mis(Trial,varargin)
% $$$ function [fet,featureTitles,featureDesc,Nmean,Nstd] = fet_mis(Trial,varargin)
% $$$ defargs = struct('newSampleRate', 12,                       ...
% $$$                  'normalize'    , false,                    ...
% $$$                  'procOpts'     , {'SPLINE_SPINE_HEAD_EQD'});



% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('newSampleRate', 12,                       ...
                 'normalize'    , false,                    ...
                 'procOpts'     , {'SPLINE_SPINE_HEAD_EQD'});

[newSampleRate,normalize,procOpts] = DefaultArgs(varargin,defargs,'--struct');
%--------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

% INIT Feature
fet = MTADfet(Trial.spath,...
              [],...
              [],...
              newSampleRate,...
              Trial.sync.copy,...
              Trial.sync.data(1),...
              [],'TimeSeries',[],'req20160310_selected_features','fet_mis','m');                  

% XYZ preprocessed 
[xyz,ss] = preproc_xyz(Trial,procOpts);
xyz.resample(newSampleRate);
ss.resample(xyz);

% XYZ filtered 
fxyz = xyz.copy;
fxyz.filter('ButFilter',3,2.4,'low');

% FVELXY Filtered marker speeds in XY plane
fvelxy = fxyz.vel({'spine_lower','spine_upper','hcom','acom'},[1,2]);
fvelxy.data(fvelxy.data<1e-4)=1e-4;
fvelxy.data = log10(fvelxy.data);

% FVELZ Filtered marker speeds in Z axis
fvelz = fxyz.vel({'hcom','acom'},[3]);
fvelz.filter('ButFilter',3,2.5,'low');
fvelz.data(fvelz.data<0)=.1;
fvelz.data = log10(fvelz.data);

% FANG Filtered Intermarker angles 
fang = create(MTADang,Trial,fxyz);

% SS 
sd = sqrt(sum((ss.data-circshift(ss.data,-1,2)).^2,3));
sn = sum(sd(:,2:end-1),2)./sd(:,end);
sv = Trial.xyz.copy;
sv.data = sn;


% AV 
sang = [circ_dist(fang(:,'spine_lower','pelvis_root',1),fang(:,'pelvis_root','spine_middle',1)),...
        circ_dist(fang(:,'pelvis_root','spine_middle',1),fang(:,'spine_middle','spine_upper',1)),...
        circ_dist(fang(:,'spine_middle','spine_upper',1),fang(:,'spine_upper','hcom',1))];
av = fang.copy;
av.data = abs(sum(sang,2)-circ_mean(sum(sang,2)));

% PPC feature
try
    man = Trial.load('fet','lsppc');
catch err
    gen_fet_lsppc(Trial);    
    man = Trial.load('fet','lsppc');
end
man.filter('ButFilter',3,2,'low');
man.resample(fxyz);

% CAT feature
fet.data = [ fang(:,'spine_lower','spine_middle',2),...
             fang(:,'spine_lower','hcom',2),        ...
             fang(:,'pelvis_root','spine_upper',2), ...
             fang(:,'spine_middle','spine_upper',2),...
             fang(:,'spine_middle','hcom',2),       ...
             man.data,                              ...
             fxyz(:,'spine_lower',3),               ...
             fxyz(:,'pelvis_root',3),               ...
             fxyz(:,'spine_middle',3),              ...             
             fxyz(:,'spine_upper',3),               ...             
             fvelxy.data,                           ...
             fvelz.data,                            ...
             sv.data,                               ...
             av.data                                ...
];

fet.data(~nniz(xyz),:)=0;
featureTitles = {};
featureDesc = {};
if nargout>1,
    % 1.
    featureTitles(end+1) = {'Pitch BLBM'};    
    %featureTitles(end+1) = {'Pitch_{BMBU}'};
    %featureTitles(end+1) = {'$\psi_{BMBU} \quad rad$'};
    featureDesc(end+1) = {['Pitch of body lower to body middle relative to xy ' ...
                           'plane']};
    %2
    featureTitles(end+1) = {'Pitch BLHC'};    
    %featureTitles(end+1) = {'Pitch_{BMBU}'};
    %featureTitles(end+1) = {'$\psi_{BMBU} \quad rad$'};
    featureDesc(end+1) = {['Pitch of body lower to head COM relative to xy ' ...
                           'plane']};
    % 3.
    featureTitles(end+1) = {'Pitch BPBU'};    
    %featureTitles(end+1) = {'Pitch_{BLBM}'};
    %featureTitles(end+1) = {'$\psi_{BLBM} \quad rad$'};
    featureDesc(end+1) = {['Pitch of body pelvis to body upper relative to xy ' ...
                           'plane']};
    %4
    featureTitles(end+1) = {'Pitch BMBU'};    
    %featureTitles(end+1) = {'Pitch_{BMBU}'};
    %featureTitles(end+1) = {'$\psi_{BMBU} \quad rad$'};
    featureDesc(end+1) = {['Pitch of body middle to body upper relative to xy ' ...
                           'plane']};
    %5
    featureTitles(end+1) = {'Pitch BMHC'};    
    %featureTitles(end+1) = {'Pitch_{BMBU}'};
    %featureTitles(end+1) = {'$\psi_{BMBU} \quad rad$'};
    featureDesc(end+1) = {['Pitch of body middle to head COM relative to xy ' ...
                           'plane']};
    % 6.
    featureTitles(end+1) = {'PPC traj yaw'};    
    %featureTitles(end+1) = {'PPC_{traj yaw}'};
    featureDesc(end+1) = {['1 Hz lowpass filtered Pair-wise Phase Consisistency(PPC) ',...
                           'of the yaw of trajectories of all makers along the ',...
                           'rostro-caudal axis']};

    % 7.
    featureTitles(end+1) = {'Z BL'};
    %featureTitles(end+1) = {'Z_{BL}'};
    featureDesc(end+1) = {['1 Hz low pass filtered height of the ' ...
                        'lower spine marker']};
    % 8
    %  '$\displaystyle\frac{dz{BL}}{dt}$'
    featureTitles(end+1) = {'Z BP'};
    %featureTitles(end+1) = {'Z_{BP}'};
    featureDesc(end+1) = {'1 Hz low pass filtered height of the pelvic maker'};

    % 9 
    featureTitles(end+1) = {'Z BM'};
    %featureTitles(end+1) = {'Z_{BM}'};
    featureDesc(end+1) = {'1 Hz low pass filtered height of the middle spine maker'};

    % 10
    featureTitles(end+1) = {'Z BU'};
    %featureTitles(end+1) = {'Z_{BU}'};
    featureDesc(end+1) = {'1 Hz low pass filtered height of the middle spine maker'};

    % 11. 
    featureTitles(end+1) = {'XY Speed BL'};
    %featureTitles(end+1) = {'d(XY_{BL})/dt'};
    %featureTitles(end+1) = {'$log10(\displaystyle\frac{dZ_{BL}}{dt}) \quad log10(mm/s)$'};
    featureDesc(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                    'the spine lower maker']};
    % 12.
    featureTitles(end+1) = {'XY Speed BU'};
    %featureTitles(end+1) = {'XY Speed_{BU}'};
    %featureTitles(end+1) = {'$log10(\displaystyle\frac{dZ_{BU}}{dt}) \quad log10(mm/s)$'};    
    featureDesc(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                    'the upper spine marker']};
    % 13.
    featureTitles(end+1) = {'XY Speed HF'};
    %featureTitles(end+1) = {'XY Speed_{HF}'};
    %featureTitles(end+1) = {'$log10(\displaystyle\frac{dZ_{HF}}{dt}) \quad log10(mm/s)$'};    
    featureDesc(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                    'the front head marker']};
    % 14.
    featureTitles(end+1) = {'XY Speed AC'};
    %featureTitles(end+1) = {'XY Speed_{AC}'};
    %featureTitles(end+1) = {'$log10(\displaystyle\frac{dZ_{AC}}{dt}) \quad log10(mm/s)$'};    
    featureDesc(end+1) = {['2.4 Hz low pass filtered speed in the xy plane of ' ...
                    'spine''s and head''s ceter of mass marker']};
    % 15.
    featureTitles(end+1) = {'Z Speed HF'};
    featureDesc(end+1) = {['1 Hz low pass filtered speed in the z axis of the ' ...
                    'front head marker']};
    % 16.
    featureTitles(end+1) = {'Z Speed AC'};
    featureDesc(end+1) = {['1 Hz low pass filtered speed in the z axis of the ' ...
                    'spine''s and head''s ceter of mass marker']};
    % 17.
    featureTitles(end+1) = {'Spine Sinuosity'};
    featureDesc(end+1) = {['Length of the spine divided by the ' ...
                        'distance beween the endpoints']};
    % 18.  
    featureTitles(end+1) = {'mean(d(yaw BLBPBMBUHC)/dt)'};
    %featureTitles(end+1) = {'d(yaw_{BLBU})/dt'};
    %featureTitles(end+1) = {'$log10(\displaystyle\frac{d\theta_{BLBU}}{dt}) \quad log10(rad/s)$'};
    featureDesc(end+1) = {['Mean Yaw speed of the vector from ' ...
                        'the lower body to the head''s center of mass']};
end


%---------------------------------------------------------------------------------------------------