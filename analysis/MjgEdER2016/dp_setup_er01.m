% Data Processing
% Project: general 
% Analysis: MjgEdER2016
% 
% subject: er01
% 
% Sessions:
%     er01-20110719
%     er01-20110721
%     er01-20110722

% LOAD / SET session meta data ------------------------------------------------------------------------
% SET parameters for each session
subject = 'er01'
clear('Sessions');
Sessions(1) = struct('sessionName',   'er01-20110719',...
                     'mazeName',      'cof',...
                     'trialName',     'all',...
                     'xyz_host',      fullfile(getenv('PROJECT'),'data','processed','xyz',subject),...
                     'nlx_host',      fullfile(getenv('PROJECT'),'data','processed','nlx'),...
                     'xyzSampleRate', 119.881035,...
                     'host',          'lmu',...
                     'TTLValue',      '0x4000'...
                     );

Sessions(end+1) = Sessions(1);
Sessions(end).sessionName = 'er01-20110721';
Sessions(end).sessionName = 'er01-20110722';
%------------------------------------------------------------------------------------------------------


% COMPILE sessions ------------------------------------------------------------------------------------
QuickSessionSetup(Sessions);
%------------------------------------------------------------------------------------------------------



% load session
s = MTASession.validate(Sessions(1));


% CORRECT rigidbody marker swaps
% ----------------------------------------------------------------------
PlotSessionErrors(s);
goodIndex = 59400;
ERCOR_fillgaps_RidgidBody(s,goodIndex,{'head_back','head_left','head_front','head_right'});
%------------------------------------------------------------------------------------------------------


% REPORT session spatial coverage----------------------------------------------------------------------
% plot timeseries for z axis
pZ(s);
% plot position within xy plane
pXY(s);









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
                    'NN0317R',...
                    'jg05-20120317.cof.all',...
                    'hl_3_jg_r',...
                    'states',{'loc','rear','pause','groom','sit'});
