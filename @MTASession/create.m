function Session = create(Session,varargin)
% Session = create(Session,varargin)
% Wrapper function for functions used to synchronize experimental
% data. The choice of function is dependent upon the combination of
% recording systems used in the session.
%  
[TTLValue,xyzSystem,ephySystem,xyzSampleRate] = DefaultArgs(varargin,{'0x8000','vicon','nlx',[]});
switch ephySystem,

  case 'nlx',
    switch xyzSystem,
      case 'vicon',
        Session = syncViconNlx(Session,TTLValue,xyzSampleRate);
    end

  case 'blackrock'
    warning(['Session creation routine does not exist ' ...
             'for Blackrock, thank you and have a nice day'])

  otherwise
    switch xyzSystem,
      case 'vicon'
        Session = loadVicon(Session,xyzSampleRate);
    end
end
end

