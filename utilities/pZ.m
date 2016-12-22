function pZ(Trial,varargin)

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
plot(xyz(:,Trial.trackingMarker,3));
try,Lines(Trial.stc{key}(:),[],'r');end