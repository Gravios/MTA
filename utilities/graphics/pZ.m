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

xyz = preproc_xyz(Trial);
figure(hfig),hold on
plot(xyz(:,'hcom',3));
try,Lines(Trial.stc{key,xyz.sampleRate}(:),[],'r');end