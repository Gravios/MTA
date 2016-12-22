function Trial = labelTheta(Trial,varargin)
%function Trial = labelTheta(Trial,varargin)
%[Stc,thetaChan,overwrite] = DefaultArgs(varargin,{[],1,false});
[Stc,thetaChan,overwrite] = DefaultArgs(varargin,{[],1,false});

if isempty(Stc),
    Stc = Trial.stc.copy;
end


sempty = isempty(Stc.gsi('t'));
if sempty||overwrite,
    if ~sempty, 
        Stc.states(Stc.gsi('t')) = [];
    end
    cpath = pwd;
    cd(Trial.spath);
    if ~exist(fullfile(Trial.spath,[Trial.name '.sts.theta']),'file')||overwrite,
        CheckEegStates(Trial.name,[],[],[],thetaChan,[],'display',overwrite);
        uiwait(gcf);
    end
end



Stc.states(Stc.gsi('t')) = [];

data = load(fullfile(Trial.spath,[Trial.name '.sts.theta']));
sync = Trial.lfp.sync.copy;
lsync = sync.sync.copy;
lsync.resample(Trial.lfp.sampleRate);
data = IntersectRanges(lsync.data,data)-lsync.data(1)+1;
Stc.addState(Trial.spath,...
             Trial.filebase,...
             data,...
             Trial.lfp.sampleRate,...
             Trial.sync.copy,...
             Trial.sync.data(1),...
             'theta','t');

Stc{'t'}.save;

Stc.save(1);
Trial.stc = Stc;
Trial.save;