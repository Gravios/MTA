function Session = create(Session,varargin)
% Session = create(Session,varargin)
% Wrapper function for functions used to synchronize experimental
% data. The choice of function is dependent upon the combination of
% recording systems used in the session.
%  
[TTLValue,dataLoggers,xyzSampleRate] = DefaultArgs(varargin,{'',{},[]});

if     all(~cellfun(@isempty,regexpi(dataLoggers,{'nlx','vicon'})))
    Session = syncViconNlx(Session,TTLValue,xyzSampleRate);

elseif all(~cellfun(@isempty,regexpi(dataLoggers,{'blackrock','vicon'})))
    warning(['Session creation routine does not exist ' ...
             'for Blackrock, thank you and have a nice day'])

elseif all(~cellfun(@isempty,regexpi(dataLoggers,{'openephys','vicon'})))
    Session = sync_openephys_vicon(Session,xyzSampleRate);
elseif all(~cellfun(@isempty,regexpi(dataLoggers,{'vicon'})))
    Session = loadVicon(Session,xyzSampleRate);

elseif all(~cellfun(@isempty,regexpi(dataLoggers,{'optitrack'})))
    Session = load_optitrack(Session,xyzSampleRate);
    
else
    warning(['MTASession:create:PatternNotFound: {', strjoin(dataLoggers,','),'}']);
end

