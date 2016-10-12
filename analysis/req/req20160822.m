
SubjectId = 'ER06';
SessionDate = '20130612';
SessionName = ['ER10-20140817']; % what was once known as filebase
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




