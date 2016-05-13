function [smat,labels,keys] = stc2mat(Stc,Data,varargin)
%function [smat,labels,keys] = stc2mat(Stc,Data)
% [states] = DefaultArgs(varargin,{Stc.list_state_attrib('key')},true);
% assumes no heirarchical relationships between states
% each is mutually exclusive from all others
if iscell(Stc), Stc = Stc{1}; end
if iscell(Data), Data = Data{1}; end

[states] = DefaultArgs(varargin,{Stc.list_state_attrib('key')},true);


Stc = Stc.copy; % Don't mess up the Stc in other workspaces

nsts = numel(states);
smat = zeros([Data.size(1),nsts]);
keys = {};
labels = {};
g = 1;

if numel(states)==1&&iscell(states)
    sts = {Stc(states{:})};    
else
    sts = Stc(states{:});        
end

for tper = sts
    tper = tper{1};
    tper.cast('TimeSeries');
    if tper.sampleRate~=Data.sampleRate,
        tper = resample(tper,Data);
    end
    smat(tper==1,g) = g;
    if nargout>1,
        keys(g) = {tper.key};
        labels(g) = {tper.label};
    end
    g = g+1;

end
end
