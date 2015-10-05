function [smat,labels,keys] = stc2mat(Stc,Data,varargin)
%function [smat,labels,keys] = stc2mat(Stc,Data)
% [states] = DefaultArgs(varargin,{Stc.list_state_attrib('key')},true);
% assumes no heirarchical relationships between states
% each is mutually exclusive from all others
[states] = DefaultArgs(varargin,{Stc.list_state_attrib('key')},true);

Stc = Stc.copy; % Don't mess up the Stc in other workspaces

nsts = numel(states);
smat = zeros([Data.size(1),nsts]);
keys = {};
labels = {};
g = 1;
for i = Stc.gsi(states), % change this for loop to something better
    tper = [Stc.states{i}];
    tper = resample(tper.cast('TimeSeries'),Data);
    smat(tper==1,g) = g;
    if nargout>1,
        keys(g) = {tper.key};
        labels(g) = {tper.label};
    end
    g = g+1;
end

end
