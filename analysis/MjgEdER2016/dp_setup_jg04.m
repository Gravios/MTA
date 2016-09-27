

S = get_session_list('jg04')

%% Setup jg05-20120311

i=1;
% setup sessions
%for i = 1:numel(S),

QuickSessionSetup(S(i));
% load session
s = MTASession.validate([S(i).sessionName,'.',S(i).mazeName,'.',S(i).trialName]);

% plot timeseries for z axis
pZ(s);
% plot position within xy plane
pXY(s);

PlotSessionErrors(s);

% correct marker swaps
ERCOR_fillgaps_RidgidBody(s,1234000);
ERCOR_fillgaps_RidgidBody(s,128800);
ERCOR_fillgaps_RidgidBody(s,328000);

PlotSessionErrors(s)

% create trial wrapper for full session
QuickTrialSetup(s)
Trial = MTATrial.validate(s.filebase);

% plot timeseries for z axisa
pZ(Trial);
% plot position within xy plane
pXY(Trial);



thetaChan = 16;
Trial = labelTheta(Trial,[],thetaChan,1);
Trial = labelTheta(Trial);

Stc = Trial.stc.copy;


overwrite = 0;
units = [];

% compute 3d place fields for the theta state
pfs_3d_theta(Trial,1,1);

