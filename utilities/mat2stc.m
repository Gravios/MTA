function [smat,labels,keys] = mat2stc(mstc,Stc,sampleRate,varargin)
%function [smat,labels,keys] = stc2mat(Stc,Data)
% assumes no heirarchical relationships between states
% each is mutually exclusive from all others
[labels,keys,sync,origin] = DefaultArgs(varargin,{Stc.list_state_attrib('key'),T},true);

Stc = Stc.copy; % Don't mess up the Stc in other workspaces

nsts = numel(states);
smat = zeros([Data.size(1),nsts]);
keys = {};
labels = {};
g = 1;
for i = Stc.gsi(states),
    tper = resample(Stc.states{i},Data);
    tper.cast('TimeSeries');
    smat(tper==1,g) = g;
    if nargout>1,
        keys(g) = {tper.key};
        labels(g) = {tper.label};
    end
    g = g+1;
end

end
