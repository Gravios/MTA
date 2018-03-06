function pXY(Trial,varargin)

switch numel(varargin)
  case 2
    hfig = varargin{1};
    if ischar(varargin{2}),
        sper = [Trial.stc{varargin{2}}];
    else,
        sper = varargin{2};
    end    
  case 1
    hfig = varargin{1};
  case 3
    hfig = varargin{1};
    if ischar(varargin{2}),
        sper = [Trial.stc{varargin{2}}];
    else,
        sper = varargin{2};
    end    
    marker = varargin{3};
    
  otherwise         
    hfig = gcf;
    sper = Trial.stc{'r'};    
    marker = Trial.trackingMarker;
end

try
    xyz = preproc_xyz(Trial,'trb');
catch
    xyz = Trial.load('xyz');
end

figure(hfig),hold on
plot(xyz(:,Trial.trackingMarker,1),xyz(:,Trial.trackingMarker,2),'.');
alims(gca,Trial.maze.visible_volume);
try,
    plot(xyz(sper,marker,1),xyz(sper,marker,2),'.r');
end