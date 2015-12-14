function printFig(Session,varargin)
%printFig(Session,varargin)
%saves an image of a target figure. 
%
%  handle     - figureHandel: default current figure handel
%
%  imFileType - string: image file type, supports 'eps' and 'png'
%
%  id         - string: identifier for the figure, (e.g. unit_number )
%
%  name       - string: name, default - random number between 
%
%  path       - string: path where the images will be saved -
%                      default ~/figures/f_"datestr"/
%
%  example:
%
%    Session.printFig(gcf,'png',num2str(1),[],[]);
%
%    Session.printFig(gcf,'png',num2str(1),'PlaceFields',Session.spath);
%
[handle,imFileType,id,name,path] = DefaultArgs(varargin,{gcf,'png',[],num2str(randi(100000,1)),['~/figures/f' datestr(now,29)]});
if ~exist(path,'dir'),
    mkdir(path),
end
if ~isempty(id), id = [ '.' num2str(id)];end
switch imFileType
  case 'eps'
    print(handle,'-dpsc2',[path '/' Session.filebase '.' name id '.' imFileType]);
  case 'png'
    print(handle,'-dpng',[path '/' Session.filebase '.' name id '.' imFileType]);
end
end
