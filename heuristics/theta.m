function [data,sampleRate,label,key] = theta(Session,varargin)
[mode,label,key] = DefaultArgs(varargin,{'sts2double','theta','t'});

switch mode
    case 'sts2double'
        sampleRate = Session.lfp.sampleRate;
        
        data = load(fullfile(Session.spath,[Session.name '.sts.theta']));
        
    case 'sts2epoch'
        data = load(fullfile(Session.spath,[Session.name '.sts.theta']));
        sync = Session.lfp.sync.copy;
        sync.resample(Session.lfp.sampleRate);
        
        data = MTADepoch(Session.spath,...
                        Session.filebase,...
                        data,...
                        Session.lfp.sampleRate,...
                        sync,...
                        sync(1),...
                        label,key);
        Session.resync(data);
        
end

end
    