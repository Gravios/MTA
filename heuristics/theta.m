function [data,sampleRate,label,key] = theta(Session,varargin)
[mode,label,key] = DefaultArgs(varargin,{'sts','theta','t'});

switch mode
    case 'sts'
        data = load(fullfile(Session.spath,[Session.name '.sts.theta']));
        sampleRate = Session.lfp.sampleRate;
end

end
    