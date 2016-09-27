
% ER09
slist = {'ER09-20131122.cof',... K F
         'ER09-20131128.cof',... K F
         'ER09-20131203.cof',... P F
         'ER09-20131204.cof'};%  P F




clear('Sessions');
Sessions(1) = struct('sessionName',   'ER09-20131122',...
                     'mazeName',      'cof',...
                     'trialName',     'all',...
                     'xyz_host',      fullfile(getenv('PROJECT'),'data','processed','xyz','ER09'),...
                     'nlx_host',      '/storage/evgeny/data/processed/ER09',...
                     'xyzSampleRate', 119.881035,...
                     'host',          'lmu',...
                     'TTLValue',      '0x0001'...
                     );

Sessions(end+1) = Sessions(1);
Sessions(end).sessionName = 'ER09-20131128';



%% Setup jg05-20120311
% setup sessions
QuickSessionSetup(Sessions);

% load session
s = MTASession('jg05-20120311');

% correct marker swaps
ERCOR_fillgaps_RidgidBody(s,1234000);

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

