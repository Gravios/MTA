function v = vel(Data,varargin)
%v = vel(Data,varargin)
%calculate the speed of marker(s)
%[markers,dims] = DefaultArgs(varargin,{[1:Data.model.N],[1:size(Data.xyz,3)]});
%zeros are added to the the beginning to for padding; reinterpolate
%later to match xyz vector.
[markers,dims] = DefaultArgs(varargin,{1:Data.model.N, 1:Data.size(3)});

if iscell(markers)||ischar(markers),
    markers = Data.model.gmi(markers);
elseif islogical(markers)
    markers = find(markers);
end

if ischar(markers),markers ={markers};end
v = MTADxyz('data',cat(1,zeros([1,numel(markers)]),sqrt(sum(diff(Data.subsref(substruct('()',{':',markers,dims})),1,1).^2,3)).*Data.sampleRate./10),...
            'sampleRate',Data.sampleRate);
v.model = Data.model.copy;
v.model.N = numel(markers);
%v.model.Connections = 
v.model.Markers = v.model.Markers(markers);
end            
