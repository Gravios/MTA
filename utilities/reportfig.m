%function reportfig(Trial, FigHandle, FileName, Preview, Comment, Resolution,SaveFig)
% This is a quick and simple function to generate an HTML
% page containing the figure and comments on it. Figures entered
% under the same filename (or on the same day) will be added to the same
% page. The html file and directory with figures are stored in 
% homedir/mrep. Both jpeg and matlab .fig are stored. 
% FigHandle  (default = gcf) 
% FileName  - (default = current date) just name, no path 
% Preview (default = 0) if 1 - launches browser to preview the page
% Comment - to put below, default - Figure # 
% SaveFig - saves the .fig file as well. default =0
function reportfig(varargin)

TodayDate = date;
DefFileName = ['mrep.' TodayDate];	
if ispc,
    DefRepPath = fullfile(getenv('HOMEPATH'),'figures');        
else
    DefRepPath = fullfile(getenv('HOME'),'figures');
end
if ~exist(DefRepPath,'dir'),mkdir(DefRepPath),end

[Trial,FigHandle, FileName, FigDir, Preview,Comment,Resolution,SaveFig,FigCount,format] = ...
    DefaultArgs(varargin,{DefRepPath,gcf,DefFileName,'',0,'',200,0,false,'png'});


%% Build Paths
myStack = dbstack;
if numel(myStack)>1,
    myFunFile = myStack(2).name;
else
    myFunFile = 'wksp';
end


if isa(Trial,'MTASession'),
    TrialFigPath= fullfile(Trial.path.data,'figures',FigDir);
    CssDir = fullfile(Trial.path.css);
elseif ischar(Trial),
    TrialFigPath = fullfile(Trial,FigDir);    
    CssDir = fullfile(fileparts(mfilename('fullpath')),'..','css');
else
    error('reportfig:BadPathError:Wait... where do you want to put it !?')
end



FunFileDir = fullfile(TrialFigPath, myFunFile);
HtmlName =   fullfile(FunFileDir,[FileName '.html']);
DirName =    fullfile(FunFileDir, FileName);



%% Create Dirs if necessary 
if ~exist(TrialFigPath,'dir'), mkdir(TrialFigPath); end
if ~exist(FunFileDir,'dir'), mkdir(FunFileDir); end


if exist(HtmlName,'file'),
    fprintf('REPORTFIG: Adding to existing report\n');
    DirCont = dir(DirName);

    DirCont([1,2]) = [];
    maxindex=1;
    for i=1:length(DirCont)
        [d1,name,ext] = fileparts(DirCont(i).name);
        digind = strfind(name,'-')+1;
        digind = digind(end);
        maxindex = max(maxindex, str2num(name(digind:end)));
    end
    FigIndex = maxindex+1;
else
    mkdir(FunFileDir,FileName);
    FigIndex = 1;
end


if isa(Trial,'MTASession'),
    if FigCount,
        FigName = [FileName '-' Trial.filebase '-' num2str(FigIndex)];
    else
        FigName = [FileName '-' Trial.filebase];
    end
else
    FigName = [FileName '-' num2str(FigIndex)];
end

LongFigName = fullfile(DirName, FigName);


%% Write HTML file 
fpHTML=fopen(HtmlName,'at+');
if fpHTML<3
    error('Can''t open html file');
end		
if FigIndex == 1
    copyfile(CssDir,fullfile(DirName,'css'));
    fprintf(fpHTML,['<html>\n<head>\n'...
                      '<link rel="stylesheet" Type="text/css" href="' FileName '/css/reportfig.css">\n',...
                      '<link rel="stylesheet" Type="text/css" href="' FileName '/css/menubar.css">',...
                      '\n<title>' ...
                      'Matlab report %s' ...
                      '</title>\n</head>\n' ...
                      '<body>\n<h2>' ...
                      'Matlab Report Page : %s</h2>\n'],...
            FileName,FileName);    
else
    fprintf(fpHTML,'<html>\n<body>\n',FileName);
end

%% Print figure 

set(FigHandle,'PaperUnits','centimeters');
set(FigHandle,'PaperPosition',[0,0,8,6]);


switch format
  case 'png'
    renderer = '-opengl';
  case 'bmp'
    renderer = '-opengl';
  case 'jpeg'
    renderer = '-opengl';
  case 'eps'
    renderer = '-painters';
  case 'pdf'
    renderer = '-painters';
end

export_fig(LongFigName,['-' format],renderer,['-r' num2str(Resolution)],FigHandle);


if SaveFig
    % and save figure
    saveas(FigHandle, [LongFigName '.fig'], 'fig');
end


% now generate link to image 
fprintf(fpHTML, '<img class="resize" src="%s" alt="Figure %d">\n<br><br>\n',[FileName '/' FigName '.' format],FigIndex);

if isempty(Comment) % Prompt for comment
    prompt = ['Enter comments for figure' num2str(FigIndex)];
    Text = inputdlg(prompt, 'Commnents dialog', 1, {['Figure ' num2str(FigIndex)]});
    Comment = Text{1};
end

% Print comment to html file
fprintf(fpHTML,'<p>%s</p>\n<br>\n<hr>\n<br>\n',Comment);


fprintf(fpHTML,'</body>\n</html>\n\n');
fclose(fpHTML);



if Preview
    web(HtmlName,'-browser');
end
