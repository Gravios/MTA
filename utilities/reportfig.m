function reportfig(varargin)
%function reportfig(varargin)
%
%  varargin:
%
%    Trial:      MTASession/string, used to specify the location where
%                figures will be saved. 
%
%                Alternate Inputs:
%                  Empty: create and use a "figures" directory in
%                         home directory
%
%                  Path:  Save figures in specified path
%
%                  MTA:   Save figures in the Session's data directory
%                                    
%    FigHandle:  graphics handle, the handle of the figure which is
%                to be reported
%
%    FileName:   string, name used to store the figure
%
%    FigDir:     string, if specified and does not exist, will create a
%                subdirectory in which the report will be loctated.
%
%    Preview:    logical, true  - open figure in browser. 
%                         false - nothing
%
%    Tag:        string, 
%
%    Comment:    string, add a comment for the caption in the web display
%
%    Resolution: integer, dpi probably 
%
%    SaveFig:    logical, true  - save a *.fig file. 
%                         false - nothing
%
%    format:     string, {'png','jpg',...} see matlab's saveas.m function
%
%    width:      numeric, final width of figure in centimeters
%
%    height:     numeric, final height of figure in centimeters
%
%    FigCount:   integer, specify the figure id (useful if you want
%                to overwrite an existing figure
%
%    InsertBreak:logical, true  - start a new line in the web display 
%                         false - nothing




%% Setup Default Arguments

TodayDate = date;
DefFileName = ['mrep.' TodayDate];	
if ispc,
    DefRepPath = fullfile(getenv('HOMEPATH'),'figures');        
else
    DefRepPath = fullfile(getenv('HOME'),'figures');
end
if ~exist(DefRepPath,'dir'),mkdir(DefRepPath),end


defArgs = {...
    ...Trial
        DefRepPath,                        ...
    ...                                    
    ...FigHandle
        gcf,                               ...
    ...
    ...FileName
        DefFileName,                       ...
    ...
    ...FigDir
        '',                                ...
    ...
    ...Preview
        false,                             ...
    ...
    ...Tag
        '',                                ...
    ...
    ...Comment
        '',                                ...
    ...
    ...Resolution
        200,                               ...
    ...
    ...SaveFig
        false,                             ...
    ...
    ...format
        'png',                             ...
    ...
    ...width
        8,                                 ...
    ...
    ...height
        6,                                 ...
    ...
    ...FigCount
        true,                              ...
    ...
    ...InsertBreak
    false
};

%% Load Default Arguments
[Trial,FigHandle, FileName, FigDir, Preview,Tag,Comment,...
 Resolution,SaveFig,format,width,height,FigCount,InsertBreak] = ...
    DefaultArgs(varargin,defArgs);

if ~iscell(format)
    format = {format};
end


%% Build Paths
myStack = dbstack;
if numel(myStack)>1,
    myFunFile = myStack(2).name;
else
    myFunFile = 'wksp';
end


if isa(Trial,'MTASession'),
    TrialFigPath= fullfile(Trial.path.data,'figures',FigDir);
    WebDir = fullfile(Trial.path.web);
elseif ischar(Trial),
    TrialFigPath = fullfile(Trial,FigDir);    
    WebDir = fullfile(fileparts(mfilename('fullpath')),'..','web');
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
    files = dir(DirName);    
    re = ['^(' FileName ')[-](\d)+[.]' format{1} '$'];% regular expression to match Session naming convention
    f = regexp({files(~cellfun(@isempty,regexp({files.name},re))).name},'[\d]+','match');
    FigIndex = max(cellfun(@str2num,cat(2,f{:})))+1;
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


%% Print figure 

set(FigHandle,'PaperUnits','centimeters');
set(FigHandle,'PaperPosition',[0,0,width,height]);


for i = 1:numel(format),
    switch format{i},
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

    export_fig(LongFigName,['-' format{i}],renderer,['-r' num2str(Resolution)],FigHandle);

end

if SaveFig
    % and save figure
    saveas(FigHandle, [LongFigName '.fig'], 'fig');
end




%% Generate/write to HTML file 

endTag = '</p></div></body></html>';

if exist(HtmlName,'file')
    htmlText = fileread(HtmlName);
    htmlText(end-numel(endTag):end) = [];
    fpHTML = fopen(HtmlName,'wt+');
    fprintf(fpHTML,htmlText);
    if InsertBreak, fprintf(fpHTML,'\n\n<br>\n\n'); end
else
    fpHTML=fopen(HtmlName,'wt+');    
end


if fpHTML<3
    error('Can''t open html file');
end		
if FigIndex == 1
    % Copy in CSS and JS files
    copyfile(WebDir,fullfile(DirName,'web'));
    fprintf(fpHTML,['<html>\n<head>\n'...
 ... JavaScript stuff
           '<script type="text/javascript" src="http://ajax.googleapis.com/ajax/libs/jquery/1.7/jquery.min.js"></script>\n',...
           '<script type="text/javascript" src="%s/web/javascript/fancybox/source/jquery.fancybox.pack.js"></script>\n\n',...
 ... CSS stuff
           '<link rel="stylesheet" href="%s/web/javascript/fancybox/source/jquery.fancybox.css"',...
                  'type="text/css" media="screen" />\n',...
           '<link rel="stylesheet" Type="text/css" href="%s/web/css/reportfig.css">\n',...
           ...'<link rel="stylesheet" Type="text/css" href="%s/web/css/menubar.css">\n\n',...
 ... Title of page
           '\n<title>' ...
           'Matlab report %s' ...
           '</title>\n\n\n</head>\n' ...
           '<body>\n<h2>' ...
           'Matlab Report Page : %s</h2>\n\n',...
 ... Initialize javascript stuff
           '<script>\n',...
                '\t$(document).ready(function() {\n',...
                    '\t\t$(''.fancybox'').fancybox();\n',...
                '\t});\n',...
           '</script>\n\n',...
 ... Container of images
           '<div class="tab_container">\n',...
           '\t<p>\n',...
                   ],...
            FileName,FileName,FileName,FileName,FileName);    
elseif InsertBreak
    fprintf(fpHTML,'\t</p>\n</div>\n\n<div class="tab_container">');
end


% Prompt for comment if none given ... no don't do that it's annoying
% $$$ if isempty(Comment) 
% $$$     prompt = ['Enter comments for figure' num2str(FigIndex)];
% $$$     Text = inputdlg(prompt, 'Commnents dialog', 1, {['Figure ' num2str(FigIndex)]});
% $$$     Comment = Text{1};
% $$$ end


% now generate link to image 
fprintf(fpHTML, '\n');
fprintf(fpHTML, '\t<div>\n');
fprintf(fpHTML, '\t\t<ul>\n');
fprintf(fpHTML, '\t\t\t<li>\n');
fprintf(fpHTML, '\t\t\t\t<a href="%s" class="fancybox" rel="gallery" title="%s">\n',[FileName '/' FigName '.' format{1}],Comment);
fprintf(fpHTML, '\t\t\t\t<img class="resize" src="%s" alt="Figure %d"/>\n\n',[FileName '/' FigName '.' format{1}],FigIndex);
fprintf(fpHTML, '\t\t\t\t</a>\n\n');
fprintf(fpHTML, '\t\t\t</li>\n\n');
% Print tag to html file
fprintf(fpHTML, '\t\t\t<li>\n');
fprintf(fpHTML, '\t\t\t\t%s\n',Tag);
fprintf(fpHTML, '\t\t\t</li>\n');
fprintf(fpHTML, '\t\t</ul>\n');
fprintf(fpHTML, '\t</div>\n');
% $$$ if InsertBreak, 
% $$$     fprintf(fpHTML,'\n\n<br>\n\n'); 
% $$$ end

fprintf(fpHTML,'%s',endTag);
fclose(fpHTML);



if Preview
    web(HtmlName,'-browser');
end
