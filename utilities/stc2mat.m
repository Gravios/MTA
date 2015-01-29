function [smat,keys,labels] = stc2mat(Stc,Data,varargin)
%function [smat,keys,labels] = stc2mat(Stc)
% assumes no heirarchical relationships between states
% each is mutually exclusive from all others
[states] = DefaultArgs(varargin,{{}},true);

if isempty(states),
    states = Stc.list_state_attrib('key');
end

nsts = numel(states);
smat = zeros([Data.size(1),nsts]);
keys = '';
labels = {};
for i = 1:nsts
    tper = resample(Stc{states{i}}.cast('TimeSeries'),Data);
    smat(tper==1,i) = i;
    if nargout>1,
        keys(i) = tper.key;
        labels{i} = tper.label;
    end
end

end
