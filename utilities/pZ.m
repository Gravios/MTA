function pZ(Trial,varargin)
if ~isempty(varargin),key = varargin{1};end
xyz = Trial.load('xyz');
figure,hold on
plot(xyz(:,Trial.trackingMarker,3));
try,Lines(Trial.stc{key}(:),[],'r');end