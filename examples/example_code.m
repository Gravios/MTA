%% How to MTA your Data
%
%% Setup        
%
% copy the following file 
%     /gpfs01/sirota/homes/share/matlab/MTA/config/MTAstartup.m
% 
%     MTAConfiguration('/gpfs01/sirota/bach/data/gravio','absolute');
%     MTAConfiguration('path_where_you_want_do_your_analysis','absolute'); 
% 
% copy the following lines to your matlab startup.m 
%     %% Labbox
%     addpath(genpath('/gpfs01/sirota/homes/share/matlab/labbox'));
%     rmpath(genpath('/gpfs01/sirota/homes/share/matlab/labbox/.git'));
%     MTAstartup
%

%% Hello Session 
% bellow is the setup info needed for playing with a few sessions

SessionName = 'Ed10-20140817'; % what was once known as filebase
xyz_path = '/gpfs01/sirota/homes/eduardo/data/xyz';
nlx_path = '/gpfs01/sirota/homes/eduardo/data/rawnlx';
MazeName = 'cof';
TrialName = 'all'; 
TTLValue = '0x0002'; % can be found in "SessionName".all.evt
overwrite = true;
xyzSampleRate = 119.881035;
ignoredViconTrials = [];
startStopShift = [0,0]; % shift boundaries of vicon starts and 
                        % stops by some number of seconds

SessionName = 'jg05-20120310';
xyz_path = '/gpfs01/sirota/data/bachdata/data/gravio/xyz';
nlx_path = '/gpfs01/sirota/data/bachdata/data/gravio/nlx';
MazeName = 'cof';
TrialName = 'all';
TTLValue = '0x0002' % can be found in "SessionName".all.evt
xyzSampleRate = 119.881035;
overwrite = true;
ignoredViconTrials = [];
startStopShift = [18,0]; % shift boundaries of vicon starts and 
                        % stops by some number of seconds




% Link data from specified path to the MTA root directory of
% the current user.
if ~isempty(xyz_path) && ~isempty(nlx_path)
    linkSession(SessionName,xyz_path,nlx_path);
end



if overwrite,

    Session = MTASession(SessionName,    ...
                         MazeName,       ...
                         overwrite,      ...
                         TTLValue,       ...
                         'xyzSampleRate', xyzSampleRate);

    
    %plots the height of the front marker on the rats head against time
    plot(Session.xyz(:,Session.trackingMarker,3));
    %plots x versus y 
    plot(Session.xyz(:,Session.trackingMarker,1),Session.xyz(:,Session.trackingMarker,2),'.');

    
    % Returns a Trial object from the data of the Session object 
    % and do automatic behavioral segmentation
    Trial = QuickTrialSetup(Session,         ...
                            TrialName,       ...
                            startStopShift,  ...
                            ignoredViconTrials);

end

%% Loading stuff

% This loads a trial
%Trial = MTATrial('Ed10-20140817','all'    ,'cof');
Trial = MTATrial (SessionName,    TrialName,MazeName);


% This performs a deep copy of an object
xyz = Trial.xyz.copy;
% This loads data into the object from file
xyz.load(Trial);


% Does the same thing as above in a single line.
xyz = Trial.load('xyz');


% You can filter data with an arbitrary window
%xyz.filter(ones([1,7])./7);
xyz.filter('ButFilter',3,50,'low'); %see help gtwin for gaussian kernals

% Load xyz without an existing Trial object... still makes a
% temporary one though so it's not super useful but it can do it
xyz = MTATrial('jg05-20120309').load('xyz').filter('ButFilter',3,50,'low');



% this is crapy though since the angles are not easily smoothed in
% the time domain, so we can create a new set based on a smoothed xyz
xyz = Trial.load('xyz');
xyz.filter('ButFilter',3,50,'low'); 
ang = create(MTADang,Trial,xyz);


% Marker segment angles are relative to the room coordinate system.
% The following code plots the distribution of head pitch during rearing
figure,hist(ang(Trial.stc{'r'},'head_back','head_front',2),100)

% Doing the same thing with numerical indexing
figure,hist(ang(Trial.stc{'r'},5,7,2),100)


%% Loading LFP

% for jg05 how to load the raw lfp of the hippocampus H64BUZ
lfp = Trial.load('lfp',1:64);

% for jg05 how to load the raw lfp of the hippocampus H32LIN
lfp = Trial.load('lfp' ,65:96);







%% Behavior Segmentation stuff - It's automatic if done by Quick trial setup
% most behaviors are stored as periods in the stc field of a Trial
% Stc should be a MTAStateCollection object.
% you can see which states are available by the following commands
disp(Trial.stc.list_state_attrib('label'))
disp(Trial.stc.list_state_attrib('key'))

% You can reference timeperiods of MTADepoch objects. These are
% stored in the Trial.stc property: stc = "State Collection".
% The periods will automatically be resampled if they don't match
% the sampling rate of the calling object.
figure,hist(ang(Trial.stc{'r'},'head_back','head_front',2),100)





%% Place Fields

% I'll comment on this another time
%units = select_units(Trial,18,'pyr');

pfs = MTAAknnpfs(Trial,...
                 'units',            [],...     % [] = all units (numbered by their clu
                 'states',           'walk',... % the state you're interested in
                 'overwrite',        true,...   % this will redo the calculations for
                                          ...     the specified units
                 'numIter',          1,...
                 'ufrShufBlockSize', 0,...
                 'binDims',          [30,30],...
                 'distThreshold',    125,...
                 'nNearestNeighbors',110);








