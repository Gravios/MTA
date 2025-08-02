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
Sessions(end+1) = Sessions(1);
Sessions(end).sessionName = 'er01-20110722';
%------------------------------------------------------------------------------------------------------


% COMPILE sessions ------------------------------------------------------------------------------------
QuickSessionSetup(Sessions);
%------------------------------------------------------------------------------------------------------


% LOAD session
s = MTASession.validate('er01-20110719.cof.all'); goodIndex = 59400;
s = MTASession.validate('er01-20110721.cof.all'); goodIndex = 312080;
s = MTASession.validate('er01-20110722.cof.all'); goodIndex = 317380;

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
