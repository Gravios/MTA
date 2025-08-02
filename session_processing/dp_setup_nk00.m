% Data Processing
% Project: general 
% Analysis: 
% 
% subject: nk00
% 
% Sessions:
%     nk00-20160315
%     

% LOAD / SET session meta data ------------------------------------------------------------------------
% SET parameters for each session
subject = 'nk00'
clear('Sessions');
Sessions(1) = struct('sessionName',   'nk00-20160315',...
                     'mazeName',      'sof',...
                     'trialName',     'all',...
                     'hostServer',   'lmu',             ...
                     'dataServer',   'lmu',             ...
                     'project',       'general',...
                     'xyz_host',      fullfile(getenv('PROJECT'),'data','processed','xyz',subject),...
                     'xyzSampleRate', 180,...
                     'host',          'lmu',...
                     'TTLValue',      'blah'...
                     );
%------------------------------------------------------------------------------------------------------


% COMPILE sessions ------------------------------------------------------------------------------------
QuickSessionSetup(Sessions);
%------------------------------------------------------------------------------------------------------


% LOAD session
s = MTASession.validate('nk00-20160315.sof.all'); goodIndex = 1700;

states = {'loc','lloc','hloc','rear','pause','lpause','hpause'};
stcMode = 'NN0317R';

% CORRECT rigidbody marker swaps
PlotSessionErrors(s);
ERCOR_fillgaps_RidgidBody(s,goodIndex,{'head_back','head_left','head_front','head_right'});
%------------------------------------------------------------------------------------------------------

% REPORT session spatial coverage----------------------------------------------------------------------
% plot timeseries for z axis
pZ(s);
% plot position within xy plane
pXY(s);

% SETUP Trial
QuickTrialSetup(s);
Trial = MTATrial(s);

% LABEL theta periods
Trial = labelTheta(Trial);

% COMPUTE 3d theta place fields
pfs_3d_theta(Trial);

% ESTIMATE neck position 
transform_rigidBody(s);

% LABLE Behavior with neural network 
labelBhv_NN(Trial,...
            stcMode,...
            'jg05-20120317.cof.all',...
            'hl_3_jg_r',...
            [],[],[],[],[],[],[],{'loc','rear','pause','groom','sit'});
Trial.load('stc',stcMode);
label_aux_bhv_reduced(Trial,stcMode,'overwrite',true);        

%hist_state_distrbutions

% COMPUTE the qualities and features of units
NeuronQuality(s,'overwrite',true);

% SELECT units for general processing
pft = pfs_2d_theta(Trial);
mrt = pft.maxRate;
units = select_units(Trial,18);
units = units(mrt(pft.data.clu(units))>1);

% COMPUTE 2d state place fields
pfs_2d_states(Trial);

% COMPUTE bootstrapped 2d state place fields
for sts = 1:numel(states),
    fprintf('process state: %s...\n',states{sts}); 
    defargs = get_default_args('MjgEdER2016','MTAAknnpfs_bs','struct');
    defargs.units = units;
    defargs.states = states{sts};
    defargs = struct2varargin(defargs);        
    pfkbs{sts} = MTAAknnpfs_bs(Trial,defargs{:});      
end










% create trial wrapper for full session
trialName = 'jg05-20120311.cof.all';
%QuickTrialSetup(trialName)
Trial = MTATrial.validate(trialName);



% plot timeseries for z axis
pZ(Trial);
% plot position within xy plane
pXY(Trial);

Trial = labelTheta(Trial,[],68,1);
Trial = labelTheta(Trial);

Stc = Trial.stc.copy;

xyz = Trial.load('xyz');
spk = Trial.spk.copy;
spk.create(Trial,xyz.sampleRate);
figure,plot(spk(10),'.')


trialName = 'jg05-20120311.cof.spk';
%QuickTrialSetup(trialName)
Trial = MTATrial.validate(trialName);
Trial = labelTheta(Trial);


% plot timeseries for z axis
pZ(Trial);
% plot position within xy plane
pXY(Trial);

overwrite = 0;
units = [];

% compute 3d place fields for the theta state
pfs_3d_theta(Trial,1,1);



labelBhv_NN(Trial,...
            stcMode,...
            'jg05-20120317.cof.all',...
            'hl_3_jg_r',...
            'states',{'loc','rear','pause','groom','sit'});


% DEMO for Laura ------------------------------------------------------------------------------------
xyz = Trial.load('xyz');
ang = create(MTADang,Trial,xyz);
fet = fet_bref(Trial,[],[],'none');
rhm = fet_rhm(Trial);
[ys,fs,ts] = fet_rhm(Trial,[],'mtchglong',0);


% Maze occupancy
figure,hold on
plot(xyz(:,'head_front',1),xyz(:,'head_front',2),'.');
xlabel('position x-axis (mm)')
ylabel('position y-axis (mm)')
title('Occupancy');

% Height of head
figure,hold on
plot([1:xyz.size(1)]./xyz.sampleRate,xyz(:,'head_front',3))
xlabel('Time (s)')
ylabel('Height (mm)')
title('Head Height During Behaivor');


% Reconstruction Error
figure,
plot([1:ang.size(1)]./ang.sampleRate,ang(:,'head_back','head_front',3))
xlabel('Time (s)');
ylabel('Inter Marker Distance (mm)');



figure();
hold('on');
plot([1:ang.size(1)]./ang.sampleRate,circ_dist(circshift(ang(:,'head_back','head_front',1),-10),...
                                               circshift(ang(:,'head_back','head_front',1),10)))
plot([1:ang.size(1)]./ang.sampleRate,circ_dist(circshift(ang(:,'spine_lower','spine_upper',1),-10),...
                                               circshift(ang(:,'spine_lower','spine_upper',1),10)))
xlabel('Time (s)');
ylabel('Angular speed (radian/sample)');
title('Angular Speed of Body and head');

% Feature Matrix
figure,
imagesc([1:xyz.size(1)]./xyz.sampleRate,1:fet.size(2),nunity(fet(:,[1:2:9,2:2:10,11:15,16:2:24,17:2:25,26:30]),[],[],[],[5,95])')
axis('xy');
caxis([-1,1]);
xlabel('Time (s)');
ylabel('Feature');
title('Body Referenced Features');

% Rhythmic Head Motion
figure();,sp = [];
sp(end+1) = subplot2(3,1,1:2,1);
imagesc(ts,fs,log10(ys.data)'),
axis('xy');
colormap('jet');
caxis([-5,-2]);
ylabel('Log10 RHM Power')
title('Sniffing Associated Rhythmic Head Motion (RHM)')
sp(end+1) = subplot2(3,1,3,1);
plot([1:rhm.size(1)]./rhm.sampleRate,rhm.data);
ylabel({'Rostro-caudal','Head Displacement (mm)'})
xlabel('Time(s)')
linkaxes(sp,'x');
