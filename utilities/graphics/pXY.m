function pXY(Trial,varargin)

switch numel(varargin)
  case 2
    hfig = varargin{1};
    key = varargin{2};
  case 1
    hfig = varargin{1};
  otherwise         
    hfig = gcf;
    key = 'r';
end

xyz = Trial.load('xyz');
figure(hfig),hold on
plot(xyz(:,Trial.trackingMarker,1),xyz(:,Trial.trackingMarker,2),'.');
alims(gca,Trial.maze.visible_volume);
try,Lines(Trial.stc{key}(:),[],'r');end