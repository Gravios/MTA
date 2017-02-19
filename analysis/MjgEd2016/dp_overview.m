
% jg05
slist = {'jg05-20120309.cof',... K
         'jg05-20120310.cof',... K 
         'jg05-20120311.cof',... K
         'jg05-20120317.cof'};%  K

% ER06
slist = {'ER06-20130611.ftm',... U
         'ER06-20130612.cof',... p
         'ER06-20130613.cof',... P
         'ER06-20130614.cof',... K F     
         'ER06-20130616.cof',... U
         'ER06-20130620.cof',... U         
         'ER06-20130624.cof',... U         
         'ER06-20130626.cof');%  U                  

% ER09
slist = {'ER09-20131122.cof',... K F
         'ER09-20131128.cof',... K F
         'ER09-20131203.cof',... P F
         'ER09-20131204.cof'};%  P F

% Ed10
slist = {'Ed10-20140812.cof',... P LD
         'Ed10-20140813.cof',... P LD         
         'Ed10-20140814.cof',... P LD         
         'Ed10-20140816.cof',... P LD                  
         'Ed10-20140817.cof',... K LDF
         'Ed10-20140818.ont',... 
         'Ed10-20140819.cof',... 
         'Ed10-20140820.cof',... 
         'Ed10-20140821.cof',... 
         'Ed10-20140822.cof',... K
         'Ed10-20140822.sov',... K F
         'Ed10-20140823.cof',...
         'Ed10-20140824.cof',... U
         'Ed10-20140825.cof'};%  U 


%% Setup jg05-20120311
% modified entry within get_session_list for jg05-20120311

% setup session
QuickSessionSetup(get_session_list('jg05-20120311'));

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

