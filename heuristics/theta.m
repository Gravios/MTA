function [data,sampleRate,label,key] = theta(Session,varargin)
[mode,label,key] = DefaultArgs(varargin,{'sts','theta','t'});

switch mode
    case 'sts'
        data = load(fullfile(Session.spath,[Session.name '.sts.theta']));
        syncPeriods = round(Session.sync([1,end]).*Session.lfp.sampleRate);
        data = IntersectRanges(data,syncPeriods)-round(Session.sync.origin*Session.lfp.sampleRate);
        data(data==0)=1;
        sampleRate = Session.lfp.sampleRate;
end

end
    