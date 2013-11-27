function varargout = MTABrowser(varargin)
% MTABROWSER MATLAB code for MTABrowser.fig
%      MTABROWSER, by itself, creates a new MTABROWSER or raises the existing
%      singleton*.
%
%      H = MTABROWSER returns the handle to a new MTABROWSER or the handle to
%      the existing singleton*.
%
%      MTABROWSER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MTABROWSER.M with the given input arguments.
%
%      MTABROWSER('Property','Value',...) creates a new MTABROWSER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before MTABrowser_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to MTABrowser_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Dependencies:
% DefaultArgs.m
% subplot2.m
%
% Packages:
% MTA 
%

% Edit the above text to modify the response to help MTABrowser

% Last Modified by GUIDE v2.5 06-May-2013 12:54:23

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @MTABrowser_OpeningFcn, ...
                   'gui_OutputFcn',  @MTABrowser_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


function MTABrowser_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to MTABrowser (see VARARGIN)

if ~ispc,javax.swing.UIManager.setLookAndFeel('com.sun.java.swing.plaf.gtk.GTKLookAndFeel');end

%% Loads Default Session 
setappdata(handles.MTABrowser,'SESSION_DATA',MTASession([]));
setappdata(handles.MTABrowser,'Session',getappdata(handles.MTABrowser,'SESSION_DATA'));
setappdata(handles.MTABrowser,'MLData',[]);

%% Initialize Sockets

%% Browser Control Sockets 

setappdata(handles.MTABrowser,'BSappSocket',{});

BSocket(handles,...
       'BSappSocket',...
        struct('control'    ,handles.BSdataManagement,...
               'app'        ,handles.DMapp,...
               'keyPressFcn',''));
BSocket(handles,...
       'BSappSocket',...
        struct('control'    ,handles.BSSetup,...
               'app'        ,handles.SUapp,...
               'keyPressFcn',''));
BSocket(handles,...
       'BSappSocket',...
        struct('control'    ,handles.BSmotionLabeling,...
               'app'        ,handles.MLapp,...
               'keyPressFcn',@(hObject,eventdata)MTABrowser('MLkeyState',hObject,eventdata,guidata(hObject))));
BSocket(handles,...
       'BSappSocket',...
        struct('control'    ,handles.BSlfpStates,...     
               'app'        ,handles.LSapp,...
               'keyPressFcn',''));

%% Motion Labeling Control Sockets
BSocket(handles,...
       'MLappSocket',...
        struct('control'    ,handles.MLdisplay,...     
               'app'        ,[],...
               'keyPressFcn','',...
               'view'       ,'MLView_display'));
BSocket(handles,...
       'MLappSocket',...
        struct('control'    ,handles.MLstates,...     
               'app'        ,handles.MLSapp,...
               'keyPressFcn','',...
               'view'       ,'MLView_states'));
BSocket(handles,...
       'MLappSocket',...
        struct('control'    ,handles.MLfeatures,...     
               'app'        ,handles.MLFapp,...
               'keyPressFcn','',...
               'view'       ,'MLView_features'));
BSocket(handles,...
       'MLappSocket',...
        struct('control'    ,handles.MLviewoptions,...     
               'app'        ,handles.MLVOapp,...  
               'keyPressFcn','',...
               'view'       ,'MLView_viewoptions')); 
       

%% Initialize Views
      
%% Motion Labeling Views                                 
BView(handles,...
      'MLapp',...
       struct('type','MLView_display',...
              'layout',{{struct('object',handles.MLviewer,...
                               'properties',struct(...
                                   'Position',[0.004768392370572207,0.005861664712778429,0.6539509536784741,0.9446952595936793],...
                                   'Visible' ,'on')),...
                        struct('object',handles.MLauxdata,...
                               'properties',struct(...
                                   'Position',[0.6614928282256478,0.005861664712778429,0.3355542095699576,0.9446952595936793],...
                                   'Visible' ,'on')),...
                        struct('object',handles.MLSapp,...
                               'properties',struct(...
                                   'Position',[0.004768392370572207,0.005861664712778429,0.3301934604904632,0.9446952595936793],...
                                   'Visible' ,'off')),...
                        struct('object',handles.MLFapp,...
                               'properties',struct(...
                                   'Position',[0.004768392370572207,0.005861664712778429,0.6532697547683922,0.9446952595936793],...
                                   'Visible' ,'off')),...
                        struct('object',handles.MLVOapp,...
                               'properties',struct(...
                                   'Position',[0.004768392370572207,0.005861664712778429,0.3301934604904632,0.9446952595936793],...
                                   'Visible' ,'off')),}})...                                            
                               );
    
BView(handles,...
      'MLapp',...
       struct('type','MLView_states',...
              'layout',{{struct('object',handles.MLviewer,...
                               'properties',struct(...
                                   'Position',[0.3399182561307902,0.005861664712778429,0.6539509536784741,0.9446952595936793],...
                                   'Visible' ,'on')),...
                        struct('object',handles.MLauxdata,...
                               'properties',struct(...
                                   'Position',[0.661492828225647812,0.005861664712778429,0.3355542095699576,0.9446952595936793],...
                                   'Visible' ,'off')),...                                
                        struct('object',handles.MLSapp,...
                               'properties',struct(...
                                   'Position',[0.004768392370572207,0.005861664712778429,0.3301934604904632,0.9446952595936793],...
                                   'Visible' ,'on')),...
                        struct('object',handles.MLFapp,...
                               'properties',struct(...
                                   'Position',[0.004768392370572207,0.005861664712778429,0.6532697547683922,0.9446952595936793],...
                                   'Visible' ,'off')),...
                        struct('object',handles.MLVOapp,...
                               'properties',struct(...
                                   'Position',[0.004768392370572207,0.005861664712778429,0.3301934604904632,0.9446952595936793],...
                                   'Visible' ,'off')),}})...                                            
                               );
    
BView(handles,...
      'MLapp',...
       struct('type','MLView_features',...
              'layout',{{struct('object',handles.MLviewer,...
                               'properties',struct(...
                                   'Position',[0.004768392370572207,0.005861664712778429,0.6539509536784741,0.9446952595936793],...
                                   'Visible' ,'off')),...
                        struct('object',handles.MLauxdata,...
                               'properties',struct(...
                                   'Position',[0.6614928282256478,0.005861664712778429,0.3355542095699576,0.9446952595936793],...
                                   'Visible' ,'on')),...                              
                        struct('object',handles.MLSapp,...
                               'properties',struct(...
                                   'Position',[0.004768392370572207,0.005861664712778429,0.3301934604904632,0.9446952595936793],...
                                   'Visible' ,'off')),...
                        struct('object',handles.MLFapp,...
                               'properties',struct(...
                                   'Position',[0.004768392370572207,0.005861664712778429,0.6532697547683922,0.9446952595936793],...
                                   'Visible' ,'on')),...
                        struct('object',handles.MLVOapp,...
                               'properties',struct(...
                                   'Position',[0.004768392370572207,0.005861664712778429,0.3301934604904632,0.9446952595936793],...
                                   'Visible' ,'off')),}})...                                            
                               );
    
BView(handles,...
      'MLapp',...
       struct('type','MLView_viewoptions',...
              'layout',{{struct('object',handles.MLviewer,...
                               'properties',struct(...
                                   'Position',[0.3399182561307902,0.005861664712778429,0.6539509536784741,0.9446952595936793],...
                                   'Visible' ,'on')),...
                        struct('object',handles.MLauxdata,...
                               'properties',struct(...
                                   'Position',[0.6614928282256478,0.005861664712778429,0.3355542095699576,0.9446952595936793],...
                                   'Visible' ,'off')),...                                 
                        struct('object',handles.MLSapp,...
                               'properties',struct(...
                                   'Position',[0.004768392370572207,0.005861664712778429,0.3301934604904632,0.9446952595936793],...
                                   'Visible' ,'off')),...
                        struct('object',handles.MLFapp,...
                               'properties',struct(...
                                   'Position',[0.004768392370572207,0.005861664712778429,0.6532697547683922,0.9446952595936793],...
                                   'Visible' ,'on')),...
                        struct('object',handles.MLVOapp,...
                               'properties',struct(...
                                   'Position',[0.004768392370572207,0.005861664712778429,0.3301934604904632,0.9446952595936793],...
                                   'Visible' ,'off')),}})...                                            
                               );


%% Load session information into selection tables
                                 
setappdata(handles.DMapp,'load_sessions',true);
guidata(hObject,handles);
setappdata(handles.MTABrowserStates,'previousBState',0);
set(handles.BSdataManagement,'Value',1);
BSdataManagement_Callback(handles.BSdataManagement, eventdata, handles)

set(handles.MLdisplay,'Value',1);
setappdata(handles.MLapp,'previousAppState',handles.MLdisplay);

% Choose default command line output for MTABrowser
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

function varargout = MTABrowser_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function BSocket(handles,socketType,appControlLink)
%% BSapp(s) are linked to control element
sockets = getappdata(handles.MTABrowser,socketType);
sockets{end+1} = appControlLink;
setappdata(handles.MTABrowser,socketType,sockets);

function BView(handles,parentApp,viewLayout)
%% BViews(s) are linked to add
views = getappdata(handles.(parentApp),'Views');
views{end+1} = viewLayout;
setappdata(handles.(parentApp),'Views',views);

function BStateChange(hObject,eventdata,handles)

if hObject~=getappdata(handles.MTABrowserStates,'previousBState');
    % State Change
    sockets = getappdata(handles.MTABrowser,'BSappSocket');
    for j = 1:length(sockets),
        if hObject==sockets{j}.control,
            setappdata(handles.MTABrowserStates,'activeSocket',sockets{j});
            set(sockets{j}.app,'Visible','on');
            set(sockets{j}.control,'Value',1);
            child = findall(sockets{j}.app);
            for k = 1:length(child),
                try
                    if child(k)~=handles.MTABrowser&&child(k)~=hObject,
                        set(child(k),'keyPressFcn',sockets{j}.keyPressFcn);
                    end
                    guidata(sockets{j}.app,child(k));
                end
            end
        else
            set(sockets{j}.app,'Visible','off');
            set(sockets{j}.control,'Value',0);
        end
    end
    setappdata(handles.MTABrowserStates,'previousBState',hObject);
else
    socket = getappdata(handles.MTABrowserStates,'activeSocket');
    set(socket.control,'Value',1);
end
guidata(hObject,handles);

function BViewChange(hObject,eventdata,handles)
activeBSocket   = getappdata(handles.MTABrowserStates,'activeSocket');
activeAppSocket = getappdata(activeBSocket.app,'activeSocket'); 
appName = get(activeBSocket.app,'Tag');

view = getappdata(activeBSocket.app,'Views');

for j = 1:length(view),
    if strcmp(activeAppSocket.view, view{j}.type),
        setappdata(activeBSocket.app,'activeView',view{j});
        for l = 1:length(view{j}.layout),
            props = fieldnames(view{j}.layout{l}.properties);
            for p = 1:length(props),
                set(view{j}.layout{l}.object,props{p},view{j}.layout{l}.properties.(props{p}));
            end
        end
    end
end

guidata(hObject,handles);

function AppStateChange(hObject,eventdata,handles)
%% App State Change
activeBSocket = getappdata(handles.MTABrowserStates,'activeSocket');
appName = get(activeBSocket.app,'Tag');

if hObject~=getappdata(handles.(appName),'previousAppState'),
    sockets = getappdata(handles.MTABrowser,[ appName 'Socket']);
    
    for j = 1:length(sockets),
        if hObject==sockets{j}.control,
            setappdata(activeBSocket.app,'activeSocket',sockets{j});
            set(sockets{j}.app,'Visible','on');
            set(sockets{j}.control,'Value',1);
            child = findall(sockets{j}.app);
            for k = 1:length(child),
                try
                    if child(k)~=handles.MTABrowser&child(k)~=hObject,
                        set(child(k),'keyPressFcn',sockets{j}.keyPressFcn);
                    end
                    guidata(sockets{j}.app,child(k));
                end
            end
            BViewChange(hObject,eventdata,handles);
            setappdata(handles.(appName),'previousAppState',hObject);
        else
            set(sockets{j}.app,'Visible','off');
            set(sockets{j}.control,'Value',0);
        end
    end
else
    socket = getappdata(activeBSocket.app,'activeSocket');
    set(socket.control,'Value',1);
end
guidata(hObject,handles);


%% DMapp code start
function BSdataManagement_Callback(hObject, eventdata, handles)
% hObject    handle to BSdataManagement (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Set Browser State
BStateChange(hObject,eventdata,handles);

%% load Sessions
if getappdata(handles.DMapp,'load_sessions'),
    Session = getappdata(handles.MTABrowser,'Session');

    % Select Valid Sessions
    files = dir(Session.path.data);
    re = '\d\d[-]\d\d\d\d\d\d\d\d$'; % regular expression
    sessionList = {files(~cellfun(@isempty,regexp({files.name},re))).name};
    
    
    
        % Recursive maze and trial search
    files = cellfun(@dir,cellfun(@strcat,repmat({Session.path.data},1,length(sessionList)), sessionList,'UniformOutput',false),'UniformOutput',false);
    mazelist = cell(1,length(files));
    triallist = cell(1,length(files));
    for i = 1:length(files),
        re = '.ses.mat$';
        mazeFileList = {files{i}(~cellfun(@isempty,regexp({files{i}.name},re))).name};
        tmazelist = {};

        for j = 1:length(mazeFileList)
            points = regexp(mazeFileList{j},'[.]');
            fileparts = [ 1 points+1; points-1 length(mazeFileList{j})]';
            if size(fileparts,1)==5,
                tmazelist{end+1} = mazeFileList{j}(fileparts(2,1):fileparts(2,2));
            end
        end

        mazelist{i} = tmazelist;
        ttriallist = {};
        for j =1:length(tmazelist),
            re = '.trl.mat$';
            trialFileList = {files{i}( ~cellfun(@isempty,...
                                               strfind({files{i}.name},tmazelist{j}))...
                                      &~cellfun(@isempty,regexp({files{i}.name},re))...
                                     ).name...
                            };
            tttriallist = {};
            for k = 1:length(trialFileList)            
                points = regexp(trialFileList{k},'[.]');
                fileparts = [ 1 points+1; points-1 length(trialFileList{k})]';
                if size(fileparts,1)==5,
                    tttriallist{end+1} = trialFileList{k}(fileparts(3,1):fileparts(3,2));
                end
            end    
            ttriallist{j} = tttriallist;
        end
        triallist{i} = ttriallist;
    end
    
    
    
   % Organize lists by subject
    re = '\d\d[-]\d\d\d\d\d\d\d\d$';
    dateList = {};
    mazeList = {};
    trialList = {};
    subjectList = {};
    subjectList{1} = sessionList{1}(1:regexp(sessionList{1},re)+1);
    tdate = {};
    tmaze = {};
    ttrial = {};
    tdate{1} = sessionList{1}(regexp(sessionList{1},re)+3:end);
    tmaze{1} = mazelist{1};
    ttrial{1} = triallist{1};
    uniqueSubject = [];

    for i = 2:length(sessionList),


        tsubject = sessionList{i}(1:regexp(sessionList{i},re)+1);
        uniqueSubject = ~cellfun(@strcmp,...
                                 subjectList,...
                                 repmat({sessionList{i}(1:regexp(sessionList{i},re)+1)},...
                                        1,length(subjectList))...
                                 );

        if length(uniqueSubject)==sum(uniqueSubject),
            if length(uniqueSubject)==sum(uniqueSubject),
                subjectList{end+1} = tsubject; %#ok<*AGROW>
            end
            dateList{end+1} = tdate;
            tdate = {};
            mazeList{end+1} = tmaze;
            tmaze = {};
            trialList{end+1} = ttrial;
            ttrial = {};

        end

        tdate{end+1} = sessionList{i}(regexp(sessionList{i},re)+3:end);
        tmaze{end+1} = mazelist{i};
        ttrial{end+1} = triallist{i};

        if i==length(sessionList),
            dateList{end+1} = tdate;
            tdate = {};
            mazeList{end+1} = tmaze;
            tmaze = {};
            trialList{end+1} = ttrial;
            ttrial = {};

        end

    end



    % Store Data in DMsessionData
    % subjectList,datelist,mazeList,trialList
    
    setappdata(handles.DMapp,'sessionList',sessionList);
    setappdata(handles.DMapp,'subjectList',subjectList);
    setappdata(handles.DMapp,'dateList',dateList);
    setappdata(handles.DMapp,'mazeList',mazeList);
    setappdata(handles.DMapp,'trialList',trialList);
    
    set(handles.SubjectList,'String',subjectList);
    set(handles.SubjectList,'Value',1);
    SubjectList_Callback(hObject, eventdata, handles);
    % ensure data is not loaded twice
    setappdata(handles.DMapp,'load_sessions',false); 

end


guidata(hObject,handles);
function DMsessionTable_CreateFcn(hObject, eventdata, handles)
function SubjectList_Callback(hObject, eventdata, handles)
% hObject    handle to SubjectList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Populate datelist based on Subject Selection
subjectIndex = get(handles.SubjectList,'Value');
dateList = getappdata(handles.DMapp,'dateList');
set(handles.DateList,'String',dateList{subjectIndex});
DateList_Callback(hObject, eventdata, handles)
function SubjectList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to SubjectList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function DateList_Callback(hObject, eventdata, handles)
% hObject    handle to DateList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

subjectIndex = get(handles.SubjectList,'Value');
dateIndex = get(handles.DateList,'Value');
mazeList = getappdata(handles.DMapp,'mazeList');
set(handles.MazeList,'String',mazeList{subjectIndex}{dateIndex});
if ~isempty(mazeList{subjectIndex}{dateIndex}),
    MazeList_Callback(hObject, eventdata, handles);
end
function DateList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to DateList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function MazeList_Callback(hObject, eventdata, handles)
% hObject    handle to MazeList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

subjectIndex = get(handles.SubjectList,'Value');
dateIndex = get(handles.DateList,'Value');
mazeIndex = get(handles.MazeList,'Value');

if mazeIndex,
     trialList = getappdata(handles.DMapp,'trialList');
%    if ~isempty(trialList{subjectIndex}{dateIndex}{mazeIndex}),
        set(handles.TrialList,...
             'String',...
             cat(2,{'full session'},trialList{subjectIndex}{dateIndex}{mazeIndex})...
            );
%     else
%         set(handles.TrialList,'String','');
%     end
end
function MazeList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MazeList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function TrialList_Callback(hObject, eventdata, handles)
% hObject    handle to TrialList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if sum(strcmp(get(handles.TrialList,'String'),''))||strcmp(get(handles.MazeList,'String'),''), 
    return;
end

subjectIndex = get(handles.SubjectList,'Value');
subjects = get(handles.SubjectList,'String');
dateIndex = get(handles.DateList,'Value');
dates= get(handles.DateList,'String');
mazeIndex = get(handles.MazeList,'Value');
mazes= get(handles.MazeList,'String');


subject = subjects{subjectIndex};
date = dates{dateIndex};
maze = mazes{mazeIndex};


% load Session
old_Session = getappdata(handles.MTABrowser,'SESSION_DATA');
if ~strcmp([subject '-' date '.' maze '.all'],old_Session.filebase),
    clear('old_Session');
    setappdata(handles.MTABrowser,'Session',MTASession([]));
    Session = MTASession([subject '-' date],maze);
    setappdata(handles.MTABrowser,'SESSION_DATA',Session);
else
    Session = getappdata(handles.MTABrowser,'SESSION_DATA');
end


trialIndex = get(handles.TrialList,'Value');
trials = get(handles.TrialList,'String');

if strcmp(trials{trialIndex},'full session'),
    setappdata(handles.MTABrowser,'Session',getappdata(handles.MTABrowser,'SESSION_DATA'));
else
    Trial = MTATrial(Session,trials{trialIndex});
    setappdata(handles.MTABrowser,'Session',Trial);
end
function TrialList_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TrialList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
%% DMapp code end


 

%% SUapp code start
function BSSetup_Callback(hObject, eventdata, handles)
% hObject    handle to BSSetup (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%% Set Browser State
BStateChange(hObject,eventdata,handles); 

if hObject ~= getappdata(handles.MTABrowserStates,'previousBState'),
    BStateChange(hObject,eventdata,handles);
end
%% SUapp code end


%% MLapp code start
function BSmotionLabeling_Callback(hObject, eventdata, handles)
% hObject    handle to BSmotionLabeling (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if hObject ~= getappdata(handles.MTABrowserStates,'previousBState'),
    
    
    MLData = getappdata(handles.MTABrowser,'MLData');
    MLData.record = {};
    MLData.recording=0;
    %set(handles.MTABrowser,'keyPressFcn','MTABrowser(''MLkeyState'',hObject,eventdata,handles)')
    Session = getappdata(handles.MTABrowser,'Session');
    xyzpos = Session.xyz.data;
    
    %% Translate connection map to pairs
    markerConnections=[];
    for i = 1:length(Session.model.Connections)
        markerConnections= cat(1,markerConnections,....
                               [Session.model.gmi(Session.model.Connections{i}.marker1),Session.model.gmi(Session.model.Connections{i}.marker2)]...
                              );
    end
    num_of_sticks= length(Session.model.Connections);    

    %% Initialization of graphical elements
    stick_options = {};
    stick_options.width = 1;
    stick_options.erase = 'normal'; %{none|xor|background}
    stick_options.line_style = '-';

    %% Plot MLxyzView
    cla(handles.MLxyzView);
    cla(handles.MLstateView);
    set(handles.MLxyzView,'XLim',[Session.maze.visible_volume(1,1) Session.maze.visible_volume(1,2)],...
                          'YLim',[Session.maze.visible_volume(2,1) Session.maze.visible_volume(2,2)],...
                          'ZLim',[Session.maze.visible_volume(3,1) Session.maze.visible_volume(3,2)],...
                          'View', [322.5 30],...
                          'Drawmode','fast',...
                          'NextPlot','add',...
                          'Projection','perspective',...
                          'XMinorGrid','on',...
                          'YMinorGrid','on');

    rotate3d(handles.MLxyzView);

    %% Lines connecting the markers
    try
        delete(MLData.sticks);
    end
    MLData.sticks = repmat({0},1,length(Session.model.Connections));
    for l=1:length(Session.model.Connections),
        MLData.sticks{l} =line('Parent',handles.MLxyzView,...
            'Tag',[Session.model.Connections{l}.marker1 ':' Session.model.Connections{l}.marker2],...
            'XData',[xyzpos(1,markerConnections(l,1),1),xyzpos(1,markerConnections(l,2),1)],...
            'YData',[xyzpos(1,markerConnections(l,1),2),xyzpos(1,markerConnections(l,2),2)],...
            'ZData',[xyzpos(1,markerConnections(l,1),3),xyzpos(1,markerConnections(l,2),3)],...
            'LineWidth', stick_options.width, ...
            ...%'Color'    , Session.model.Connections{l}.color./255, ...
            'EraseMode', stick_options.erase,...
            'Visible','on');
        guidata(MLData.sticks{l},handles);
    end
    
    
    %% Markers
    marker_options = {};
    marker_options.style ='o' ;
    marker_options.size =6 ;
    marker_options.erase = 'normal'; %{none|xor|background}
    try
        delete(MLData.markers);
    end
    MLData.markers = repmat({0},1,Session.model.N);
    for l=1:Session.model.N,
        MLData.markers{l} =line('Parent',handles.MLxyzView, ...
            'Tag',Session.model.Markers{l}.name,...
            'XData',[xyzpos(1,l,1),xyzpos(1,l,1)],...
            'YData',[xyzpos(1,l,2),xyzpos(1,l,2)],...
            'ZData',[xyzpos(1,l,3),xyzpos(1,l,3)],...
            'Marker'         , marker_options.style, ...
            'MarkerEdgeColor', Session.model.Markers{l}.color./255, ...
            'MarkerSize'     , marker_options.size, ...
            'MarkerFaceColor', Session.model.Markers{l}.color./255, ...
            'EraseMode'      , marker_options.erase,...
            'Visible','on');
        guidata(MLData.markers{l},handles);
    end
    
    
    %% Plot MLstateView
    %% Set labels & keys
    % Convert Bhv.States.state from Periods to time-series
    States = Session.stc.states;
    if ~isempty(States)
        labels = {};
        keys = '';
        for i = 1:length(States)
            labels{i} = States{i}.label;
            keys = strcat(keys,States{i}.key);
            tmpState = States{i}.data;
            tmpState(tmpState==0) = 1;
            States{i}.data = zeros(Session.xyz.size(1),1);
            for j = 1:size(tmpState,1),
                States{i}.data(tmpState(j,1):tmpState(j,2)) = 1;
            end
        end
    else
        labels = {'default'};
        keys = 'd';
        States{1}.data = zeros(Session.xyz.size(1),1);
    end

    current_label = 1;
    num_of_states = length(States);
    selected_states = true(1,num_of_states);
    f_line_colors = jet(num_of_states);
    
    % States,num_of_states,sxyz

    %% Set Up States MLstateView Window
    sxyz = size(xyzpos);
    for l=1:num_of_states
        States{l}.data(States{l}.data(:)>0) =States{l}.data(States{l}.data(:)>0)-1+l/2;
    end
    statesRange = [-0.5,num_of_states/2+0.5];
    windowSize = 100; 
    index_range = zeros(sxyz(1),2);
    index_range(1:windowSize,1) = 1;
    index_range(1:windowSize,2) = (1+windowSize):2*windowSize;
    index_range(windowSize+1:sxyz(1)-windowSize,1)=1:(sxyz(1)-2*windowSize);
    index_range(windowSize+1:sxyz(1)-windowSize,2)=(1+2*windowSize):(sxyz(1));
    index_range(sxyz(1)-windowSize+1:sxyz(1),1)=(sxyz(1)-2*windowSize+1):sxyz(1)-windowSize;
    index_range(sxyz(1)-windowSize+1:sxyz(1),2)=sxyz(1);

    %% Draws the Behavioral States
    state_line_options= {};
    state_line_options.width = 1;
    state_line_options.erase = 'normal'; %{none|xor|background}
    state_line_options.line_style = '-';

   
    for l=1:num_of_states,
        MLData.state_lines{l} =line('Parent', handles.MLstateView,...
                          'Tag',       labels{l},...
                          'XData',     1:windowSize,...
                          'YData',     States{l}.data(1:windowSize),...
                          'LineWidth', state_line_options.width, ...
                          'Color'    , f_line_colors(l,:), ...
                          'EraseMode', state_line_options.erase,...
                          'LineStyle', state_line_options.line_style,...
                          'Visible',   'on');
        guidata(MLData.state_lines{l},handles);
    end
    
    
    
    set(handles.MLstateView,'Xlim',[-windowSize,windowSize],...
                            'YLim',[statesRange(1),statesRange(2)],...
                            'YTick',0.5:0.5:num_of_states/2,...
                            'YTicklabel',strcat(labels, ...
                                                repmat({' ('},1,length(labels)),...
                                                str2cell(keys)',...
                                                repmat({')'},1,length(labels))) ,...
                            'YTickMode','manual');
    
    MLData.state_centerline = line([1,1], [statesRange(1),statesRange(2)],...
                           'Tag','state_centerline',...
                           'Color','r',...
                           'LineStyle','--',...
                           'LineWidth',1,...
                           'Parent', handles.MLstateView);
    guidata(MLData.state_centerline,handles)

    %% Initialize MLauxData

    try      
        delete(MLData.feature_display_handles);
        delete(MLData.feature_display_axes);
        delete(MLData.feature_display_centerlines);
        delete(MLData.feature_display_YLimCB);       
        delete(MLData.feature_display_slider);       
        delete(MLData.feature_display_panels);              
    end
    
    if ~isempty(Session.fet),
        MLData.Features = Session.fet.Features;
    else
        MLData.Features = {};
    end
    MLData.num_of_features = length(MLData.Features);
    if ~isempty(MLData.Features)
        MLData.featurelabels = {};
        for i = 1:MLData.num_of_features
            MLData.featurelabels{i} = MLData.Features{i}.label;            
        end
    else
        MLData.featurelabels = {'default'};
        MLData.Features{1} = MTAFeature('default','zeros(Session.xyz.size(1),1);',0,'line');
        MLData.num_of_features = 1;
    end
    
    
    MLData.selected_features = false(1,MLData.num_of_features);
    
    feature_line_options= {};
    feature_line_options.width = 1;
    feature_line_options.erase = 'normal'; %{none|xor|background}
    feature_line_options.line_style = '-';
    MLData.feature_line_options = feature_line_options;
    
    MLData.feature_display_panels =      repmat({0},1,MLData.num_of_features);
    MLData.feature_display_YLimCB =      repmat({0},1,MLData.num_of_features);
    MLData.feature_display_handles =     repmat({0},1,MLData.num_of_features);
    MLData.feature_display_axes =        repmat({0},1,MLData.num_of_features);
    MLData.feature_display_slider =      repmat({0},1,MLData.num_of_features);
    MLData.feature_display_centerlines = repmat({0},1,MLData.num_of_features);
    %% Animation Init Vars
    
    MLData.idx=1;
    MLData.flag_tag = false;
    MLData.erase_tag = false;
    MLData.animation = false;
    MLData.play_speed = 1;

    MLData.line_options = stick_options;
    
    MLData.xyzpos = xyzpos;
    MLData.sessionLength = Session.xyz.size(1);
    MLData.markerConnections = markerConnections;
    MLData.num_of_sticks = num_of_sticks;
    MLData.num_of_markers = Session.model.N;

    
    % MLstateView Stuff
    
    MLData.States = States;
    MLData.labels = labels;
    MLData.keys = keys;
    MLData.windowSize = 100;
    MLData.current_label = current_label;
    MLData.previous_label = labels{current_label};

    MLData.state_line_options = state_line_options;
    MLData.statesRange = statesRange;
    MLData.num_of_states = num_of_states;
    MLData.selected_states = selected_states;
    MLData.index_range = index_range;

    
    MLData.t = d2t(Session.xyz.size(1),Session.xyz.sampleRate,0);
    
    set(handles.MLposition_slider,'Value',1,'Min',1,'Max',MLData.sessionLength-MLData.windowSize, 'SliderStep', [10 100]./MLData.sessionLength);

    setappdata(handles.MLapp,'play_speed',1);
    setappdata(handles.MLapp,'animation',0);
    setappdata(handles.MLapp,'paused',0);
    setappdata(handles.MTABrowser,'MLData',MLData);
    
    %% Set Browser State
    BStateChange(hObject,eventdata,handles);     
    
    setappdata(handles.MLapp,'previousAppState',0);
    
    if hObject ~= getappdata(handles.MLapp,'previousAppState'),
        AppStateChange(handles.MLdisplay,eventdata,handles);
    end
    
    guidata(hObject, handles);
    
else
    %% Set Browser State
    BStateChange(hObject,eventdata,handles);
    guidata(hObject, handles);
end


%% MLxyzView Controls - start

% --- Executes on button press in MLviewRecord.
function MLviewRecord_Callback(hObject, eventdata, handles)
selected = get(hObject,'Value');
MLData = getappdata(handles.MTABrowser,'MLData');    
if selected,    
    MLData.record = {};
    MLData.recording=1;
else
    Session = getappdata(handles.MTABrowser,'Session');
    MLData.recording=0;
    record = MLData.record;
    save(fullfile(Session.spath, [Session.filebase ...
        '.frameset_' ...
        num2str(record{1}.index) '_' ...
        num2str(record{end}.index) '.mat']),'record');
end
setappdata(handles.MTABrowser,'MLData',MLData);

function MLposition_slider_Callback(hObject, eventdata, handles)
% hObject    handle to MLposition_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
idx = round(get(hObject,'Value'));
MLData = getappdata(handles.MTABrowser,'MLData');
MLData.idx = idx;
updateCircle(hObject,handles,0)
setappdata(handles.MTABrowser,'MLData',MLData);
set(hObject,'Value',MLData.idx);
guidata(hObject,handles);
function MLplayspeed_slider_Callback(hObject, eventdata, handles)
% hObject    handle to MLplayspeed_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.MLapp,'play_speed',get(hObject,'Value'));
set(hObject,'Value', getappdata(handles.MLapp,'play_speed'));
function MLplayforward_Callback(hObject, eventdata, handles)
% hObject    handle to MLplayforward (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.MLapp,'animation',~getappdata(handles.MLapp,'animation'));
setappdata(handles.MLapp,'paused',0);
while getappdata(handles.MLapp,'animation')
    updateCircle(hObject,handles,1);
end
function MLjumpforward_Callback(hObject, eventdata, handles)
MLData = getappdata(handles.MTABrowser,'MLData');
current_label = MLData.current_label;
setappdata(handles.MLapp,'animation',0);
setappdata(handles.MLapp,'paused',0);
if current_label~=0
    current_state = MLData.States{current_label}.data(MLData.idx);
    idx_shift = find(MLData.States{current_label}.data(MLData.idx:end)~=current_state,1);
    if ~isempty(idx_shift),
        MLData.idx = MLData.idx + idx_shift;
        MLData.flag_tag = 0;
        erase_tag = MLData.erase_tag;
        MLData.erase_tag = 0;
        setappdata(handles.MTABrowser,'MLData',MLData),
    end
    updateCircle(hObject,handles,0)
    MLData = getappdata(handles.MTABrowser,'MLData');
    MLData.flag_tag = 1;
    MLData.erase_tag = erase_tag;
    setappdata(handles.MTABrowser,'MLData',MLData),
else
    updateCircle(hObject,handles,50)
end
function MLplayreverse_Callback(hObject, eventdata, handles)
% hObject    handle to MLplayreverse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
setappdata(handles.MLapp,'animation',~getappdata(handles.MLapp,'animation'));
setappdata(handles.MLapp,'paused',0);
while getappdata(handles.MLapp,'animation')
    updateCircle(hObject, handles,-1);
end
function MLjumpback_Callback(hObject, eventdata, handles)
MLData = getappdata(handles.MTABrowser,'MLData');
current_label = MLData.current_label;
setappdata(handles.MLapp,'animation',0);
setappdata(handles.MLapp,'paused',0);
if current_label~=0
    current_state = MLData.States{current_label}.data(MLData.idx);
    idx_shift = find(MLData.States{current_label}.data(1:MLData.idx)~=current_state,1,'last');
    if ~isempty(idx_shift),
        MLData.idx = MLData.idx + idx_shift - MLData.idx;
        MLData.flag_tag = 0;
        erase_tag = MLData.erase_tag;
        MLData.erase_tag = 0;
        setappdata(handles.MTABrowser,'MLData',MLData),
    end
    updateCircle(hObject,handles,0)
    MLData = getappdata(handles.MTABrowser,'MLData');
    MLData.flag_tag = 1;
    MLData.erase_tag = erase_tag;
    setappdata(handles.MTABrowser,'MLData',MLData),
else
    updateCircle(hObject,handles,-50)
end
function MLindexinput_Callback(hObject, eventdata, handles)
% hObject    handle to MLindexinput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MLData = getappdata(handles.MTABrowser,'MLData');
MLData.idx = str2double(get(hObject,'String'));
setappdata(handles.MTABrowser,'MLData',MLData);
function MLtimeinput_Callback(hObject, eventdata, handles)
% hObject    handle to MLtimeinput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MLData = getappdata(handles.MTABrowser,'MLData');
MLData.idx = NearestNeighbor(MLData.t,str2double(get(hObject,'String')));
setappdata(handles.MTABrowser,'MLData',MLData);

%CreateFcn's
function MLtimeinput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MLtimeinput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function MLposition_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MLposition_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject, 'Value', 1);
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function MLindexinput_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MLindexinput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function MLplayspeed_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MLplayspeed_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject, 'Value', 32);
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
function MLjumpback_CreateFcn(hObject, eventdata, handles)
icon = imread('skip_backward.png');
set(hObject,'CData',icon);
guidata(hObject, handles);
function MLplayreverse_CreateFcn(hObject, eventdata, handles)
icon = imread('reverse.png');
set(hObject,'CData',icon);
guidata(hObject, handles);
function MLplayforward_CreateFcn(hObject, eventdata, handles)
icon = imread('play.png');
set(hObject,'CData',icon);
guidata(hObject, handles);
function MLjumpforward_CreateFcn(hObject, eventdata, handles)
icon = imread('skip_forward.png');
set(hObject,'CData',icon);
guidata(hObject, handles);
%% end - MLxyzView Controls



%% MLapp Display 
function MLdisplay_Callback(hObject, eventdata, handles)
AppStateChange(hObject,eventdata,handles);



%% MLS state Functions start
function MLstates_Callback(hObject, eventdata, handles)
AppStateChange(hObject,eventdata,handles);
MLSstates_loadTable(hObject, eventdata, handles);
function MLstateView_refresh(hObject,eventdata,handles)
%% Set Up States MLstateView Window
MLData = getappdata(handles.MTABrowser,'MLData');

selected_states = find(MLData.selected_states);
num_of_selected_states = length(selected_states);
f_line_colors = jet(num_of_selected_states);
f_line_colors = (f_line_colors+repmat([0.1,0,0],size(f_line_colors,1),1))./1.1;
%% Draws the Behavioral States
ssl = 1;
for l=selected_states,
    MLData.States{l}.data(MLData.States{l}.data(:)>0)=MLData.States{l}.data(MLData.States{l}.data(:)>0)./max(MLData.States{l}.data)-1+ssl/2;
    if l>length(MLData.state_lines),
        MLData.state_lines{l} =line('Parent', handles.MLstateView,...
            'XData',     1:MLData.windowSize+1,...
            'YData',     MLData.States{l}.data(MLData.idx:MLData.idx+MLData.windowSize),...
            'LineWidth', MLData.state_line_options.width, ...
            'Color'    , f_line_colors(ssl,:), ...
            'EraseMode', MLData.state_line_options.erase,...
            'LineStyle', MLData.state_line_options.line_style,...
            'Visible',   'on');   
        guidata(MLData.state_lines{l},handles)
    else
        set(MLData.state_lines{l},'Visible','on');
    end
    ssl = ssl+1;
end
unselected_states = setdiff(1:length(MLData.state_lines),selected_states);
set(cell2mat(MLData.state_lines(unselected_states)),'Visible','off');

MLData.statesRange = [-0.5,num_of_selected_states/2+0.5];
set(handles.MLstateView,'Xlim',[-MLData.windowSize,MLData.windowSize],...
    'YLim',[MLData.statesRange(1),MLData.statesRange(2)],...
    'YTick',0.5:0.5:num_of_selected_states/2,...
    'YTicklabel',strcat(MLData.labels(selected_states), ...
                        repmat({' ('},1,length(selected_states)),...
                        str2cell(MLData.keys(selected_states))',...
                        repmat({')'},1,length(selected_states))) ,...
    'YTickMode','manual');

set(MLData.state_centerline,'YData',[MLData.statesRange(1),MLData.statesRange(2)]);


if find(selected_states==MLData.current_label),
    for l = selected_states,
        if strcmp(MLData.labels{l},MLData.previous_label)
            MLData.current_label = find(selected_states==l);
        end
    end
    if isempty(MLData.current_label)
        MLData.current_label = 0;
    end
else
    MLData.current_label = 0;
end
setappdata(handles.MTABrowser,'MLData',MLData);

guidata(hObject, handles);
updateCircle(hObject,handles,0)

function MLSstateTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to MLSstateTable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
MLData = getappdata(handles.MTABrowser,'MLData');
col = eventdata.Indices(2);
row = eventdata.Indices(1);
stateTableData = get(hObject,'Data');
switch col
    case 1
        if isempty(eventdata.NewData)|sum(~cellfun(@isempty,regexpi(MLData.labels,eventdata.NewData))),
            warndlg({'New state name field is empty','or the name is arleady set','for another state'}),
            stateTableData{row,col} = eventdata.PreviousData;
            set(hObject,'Data',stateTableData);
            return
        end
        MLData.labels{row} = eventdata.NewData;
        set(handles.MLSnewStateName,'String',eventdata.NewData);
    case 2        
        if isempty(eventdata.NewData)|~isempty(find(MLData.keys==eventdata.NewData,1)),
            warndlg({'New state key field is empty','or the key is arleady set','for another state'}),
            stateTableData{row,col} = eventdata.PreviousData;
            set(hObject,'Data',stateTableData);
            return
        end
        MLData.keys(row) = eventdata.NewData;
        set(handles.MLSnewStateKey,'String',eventdata.NewData);
    case 3
        MLData.selected_states(row) = ~MLData.selected_states(row);
    otherwise
        return
end
setappdata(handles.MTABrowser,'MLData',MLData);
guidata(hObject, handles);
MLstateView_refresh(hObject,eventdata,handles);
function MLSstateTable_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to MLSstateTable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(eventdata.Indices);
    row = eventdata.Indices(1);
    stateTableData = get(handles.MLSstateTable,'Data');
    set(handles.MLSnewStateName,'String',stateTableData{row,1});
    set(handles.MLSnewStateKey,'String',stateTableData{row,2});
    guidata(hObject, handles);
end
function MLSstates_loadTable(hObject, eventdata, handles)
MLData = getappdata(handles.MTABrowser,'MLData');
columnname =   {'State','keybinding','visible','total count'};
columnformat = {'char' ,'char'      ,'logical','numeric'};
columneditable=[true  , true      , true    ,false];
statesTableData = {};
for i = 1:length(MLData.labels),
    statesTableData{i,1} = MLData.labels{i};
    statesTableData{i,2} = MLData.keys(i);
    statesTableData{i,3} = MLData.selected_states(i);
    statesTableData{i,4} = size(ThreshCross(MLData.States{i}.data,0.0005,1),1);
end

set(handles.MLSstateTable,'ColumnName',columnname);
set(handles.MLSstateTable,'ColumnFormat',columnformat);
set(handles.MLSstateTable,'ColumnEditable',columneditable);
set(handles.MLSstateTable,'RowName',[]);
set(handles.MLSstateTable,'Data',statesTableData);
set(handles.MLSstateTable,'ColumnWidth',{147,80,80,150})
guidata(handles.MLSstateTable, handles);

function MLSaddState_Callback(hObject, eventdata, handles)
% Get the name and key of the state which will be added
MLData = getappdata(handles.MTABrowser,'MLData');
newStateName = get(handles.MLSnewStateName,'String');
if isempty(newStateName)|sum(~cellfun(@isempty,regexpi(MLData.labels,newStateName))),
    warndlg({'New state name field is empty','or the name is arleady set','for another state'}),
    return
end
newStateKey = get(handles.MLSnewStateKey,'String');
if isempty(newStateKey)|~isempty(find(MLData.keys==newStateKey,1)),
    warndlg({'New state key field is empty','or the key is arleady set','for another state'}),
    return
end


MLData.labels{end+1} = newStateName;
MLData.keys(end+1) = newStateKey;
MLData.States{end+1} = MTAState(newStateKey,newStateName,zeros(length(MLData.xyzpos),1));
MLData.num_of_states = MLData.num_of_states + 1;
MLData.selected_states(end+1) = true;
setappdata(handles.MTABrowser,'MLData',MLData);
guidata(hObject, handles);
MLSstates_loadTable(hObject, eventdata, handles);
MLstateView_refresh(hObject,eventdata,handles);
function MLSremoveState_Callback(hObject, eventdata, handles)
MLData = getappdata(handles.MTABrowser,'MLData');
newStateName = get(handles.MLSnewStateName,'String');
removeStateIndex = find(~cellfun(@isempty,regexpi(MLData.labels,newStateName)));
keptStates = setdiff(1:MLData.num_of_states,removeStateIndex);
MLData.labels = MLData.labels(keptStates);
MLData.keys = MLData.keys(keptStates);
MLData.States = MLData.States(keptStates);
MLData.num_of_states = MLData.num_of_states - 1;
MLData.selected_states = MLData.selected_states(keptStates);
delete(MLData.state_lines{removeStateIndex});
MLData.state_lines = MLData.state_lines(keptStates);
setappdata(handles.MTABrowser,'MLData',MLData);
guidata(hObject, handles);
MLSstates_loadTable(hObject, eventdata, handles);
MLstateView_refresh(hObject,eventdata,handles);
function MLSmoveStateUp_Callback(hObject, eventdata, handles)
function MLSmoveStateDown_Callback(hObject, eventdata, handles)
function MLSbhvOpen_Callback(hObject, eventdata, handles)
Session = getappdata(handles.MTABrowser,'Session');
MLData = getappdata(handles.MTABrowser,'MLData');
if isempty(Session.Bhv),
    bhv_name = 'default';
else
    bhv_name = Session.Bhv.mode;
end

[filename,filepath] = uigetfile(['*.' Session.trialName '.bhv.*'],'Open Behavior Collection:',...
                     fullfile(Session.spath, [Session.filebase '.bhv.' bhv_name '.mat']),...
                     'MultiSelect','off');
                 
if ~ischar(filename)||~ischar(filepath),                 
    return
else                 
    load([filepath filename]);
    Session.Bhv = Bhv;
    States = Bhv.States;
    if ~isempty(States)
        labels = cell(1,length(States));
        keys = '';
        for i = 1:length(States)
            labels{i} = States{i}.label;
            keys = strcat(keys,States{i}.key);
            tmpState = States{i}.data;
            tmpState(tmpState==0) = 1;
            States{i}.data = zeros(Session.xyz.size(1),1);
            for j = 1:size(tmpState,1),
                States{i}.data(tmpState(j,1):tmpState(j,2)) = 1;
            end
        end        
    end
    MLData.States = States;
    MLData.labels = labels;
    MLData.keys = keys;
    MLData.current_label = 1;
    MLData.previous_label = labels{MLData.current_label};
    MLData.num_of_states = length(States);
    MLData.selected_states = true(1,MLData.num_of_states);
    setappdata(handles.MTABrowser,'MLData',MLData)
    set(handles.MLStitle,'String',['Bhavioral States - ' filename]);
    guidata(hObject, handles);
    MLSstates_loadTable(hObject, eventdata, handles);
    MLstateView_refresh(hObject,eventdata,handles);
end
function MLSbhvSave_Callback(hObject, eventdata, handles)
Session = getappdata(handles.MTABrowser,'Session');
MLData = getappdata(handles.MTABrowser,'MLData');
if isempty(Session.Bhv),
    bhv_name = 'default';
else
    bhv_name = Session.Bhv.mode;
end
[filename,filepath] = uiputfile(['*' Session.trialName '.bhv.*'],'Save Behavior Collection:',...
                     fullfile(Session.spath, [Session.filebase '.bhv.' bhv_name '.mat']));

if ~ischar(filename)||~ischar(filepath),
    return
else
    points = regexp(filename,'[.]');
    fileparts = [1 points+1; points-1 length(filename)]';
    if size(fileparts,1)==6,
        bhv_name = filename(fileparts(5,1):fileparts(5,2));
    else
        error('Bhv filename does not match standard format,')
    end
    if ~exist([filepath filename],'file'),
        Bhv = MTABhv(Session,bhv_name,1);
    else
        load([filepath filename]);
    end
    
    Bhv.States = {};
    for s = 1:MLData.num_of_states,
        state = MLData.States{s}.data./max(MLData.States{s}.data);
        state = ThreshCross(state,0.5,1);
        Bhv = Bhv.addState(MLData.keys(s),MLData.labels{s},state);
    end
    Bhv.save(Session,1);
    Session.Bhv = Bhv;
    Session.save();
end


function MLSaddState_CreateFcn(hObject, eventdata, handles)
icon = imread('add.png');
set(hObject,'Cdata',icon);
guidata(hObject, handles);
function MLSremoveState_CreateFcn(hObject, eventdata, handles)
icon = imread('subtract.png');
set(hObject,'Cdata',icon);
guidata(hObject, handles);
function MLSnewStateName_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function MLSnewStateKey_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function MLSmoveStateDown_CreateFcn(hObject, eventdata, handles)
icon = imread('arrow_down.png');
set(hObject,'CData',icon);
guidata(hObject, handles);
function MLSmoveStateUp_CreateFcn(hObject, eventdata, handles)
icon = imread('arrow_up.png');
set(hObject,'Cdata',icon);
guidata(hObject, handles);
function MLSbhvOpen_CreateFcn(hObject, eventdata, handles)
icon = imread('open.png');
set(hObject,'Cdata',icon);
guidata(hObject, handles);
function MLSbhvSave_CreateFcn(hObject, eventdata, handles)
icon = imread('save.png');
set(hObject,'Cdata',icon);
guidata(hObject, handles);

%% MLS state Functions end



%% MLF feature Functions start
function MLfeatures_Callback(hObject, eventdata, handles)
AppStateChange(hObject,eventdata,handles);
MLFfeatures_loadTable(hObject, eventdata, handles);
function MLfeatureView_refresh(hObject,eventdata,handles)
%% Set Up States MLstateView Window
MLData = getappdata(handles.MTABrowser,'MLData');

selected_features = find(MLData.selected_features);
num_of_selected_features = length(selected_features);

display_panel_position = ...
    [0.002*ones(num_of_selected_features,1),...
    1-(1/num_of_selected_features:1/num_of_selected_features:1)'+0.018405,...
    0.9950204498977505*ones(num_of_selected_features,1),...
    (1/(num_of_selected_features+.2))*ones(num_of_selected_features,1)];

display_axes_position = ...
    [0.050100502512562814,...
     0.07009345794392521,...
     0.8677219430485762,...
     0.8124299065420556];

display_slider_position = ...
    [0.9246231155778895,...
     0.07009345794392521,...
     0.05360134003350084,...
     0.8224299065420556];
 
 display_YLimCB_position = ...
    [0.9246231155778895,...
     0.90009345794392521,...
     0.04224134003350084,...
     0.08224299065420556];
 

%% Draws the Behavioral States
ssl = 1;
for l=selected_features,
    if isempty(MLData.Features{l}.feature)|...
       l>length(MLData.feature_display_handles)|...
       MLData.feature_display_axes{l}==0|...
       MLData.feature_display_handles{l}==0|...
       MLData.num_of_features==1,
 
        if isempty(MLData.Features{l}.feature),
            Session = getappdata(handles.MTABrowser,'Session');
            MLData.Features{l} = MLData.Features{l}.compute(Session);
        end
        
        if l>length(MLData.feature_display_handles)|MLData.feature_display_panels{l}==0,
            MLData.feature_display_panels{l} = uipanel('Parent',handles.MLauxdata,...
                'Title',MLData.featurelabels{l},...
                'Units','normalized',...    
                'Position',display_panel_position(ssl,:));
            
            guidata(MLData.feature_display_panels{l},handles);
            
            MLData.feature_display_axes{l}=axes('Position',display_axes_position,...
                'Units','normalized',...        
                'Parent', MLData.feature_display_panels{l},...                
                'XLim',   [-MLData.windowSize,MLData.windowSize],...
                'YLimMode','auto',...
                'FontSize',6,...
                'DrawMode',  'fast');
            guidata(MLData.feature_display_axes{l},handles);            
            
            MLData.feature_display_centerlines{l} = line([1,1], ...
                ylim(MLData.feature_display_axes{l}),...
                'Color','r',...
                'LineStyle','--',...
                'LineWidth',1,...
                'Parent', MLData.feature_display_axes{l});
            guidata(MLData.feature_display_centerlines{l},handles);
            
            MLData.feature_display_slider{l} = uicontrol(MLData.feature_display_panels{l},...
                'Units','normalized',...
                'Tag','feature_display_slider',...
                'Position',display_slider_position,...
                'Style','slider',...
                'SliderStep',[0.01,0.1],...
                'Visible','off',...
                'Callback',@(hObject,eventdata)MTABrowser('MLauxdataDisplayOption_Callback',hObject,eventdata,guidata(hObject)),...
                'Min',min(MLData.Features{l}.feature(:)),...                
                'Max',max(MLData.Features{l}.feature(:))+1,...
                'Value',mean(MLData.Features{l}.feature(:))+1);                 
            guidata(MLData.feature_display_slider{l},handles);
            
            MLData.feature_display_YLimCB{l} = uicontrol(MLData.feature_display_panels{l},...
                'Style','checkbox',...
                'Units','normalized',...        
                'Position',display_YLimCB_position,...
                'Tag', 'feature_display_YLimCB',...
                'Callback',@(hObject,eventdata)MTABrowser('MLauxdataDisplayOption_Callback',hObject,eventdata,guidata(hObject)),...
                'Visible',   'off');
            guidata(MLData.feature_display_YLimCB{l},handles);
        end
        
        switch MLData.Features{l}.plot_type,
            case 'line'
                if size(MLData.Features{l}.feature,2)==1,
                    MLData.feature_display_handles{l} =line('Parent', MLData.feature_display_axes{l},...
                        'XData',     1:MLData.windowSize+1,...
                        'YData',     MLData.Features{l}.feature(MLData.idx:MLData.idx+MLData.windowSize),...
                        'LineWidth', MLData.state_line_options.width, ...
                        'EraseMode', MLData.state_line_options.erase,...
                        'LineStyle', MLData.state_line_options.line_style,...
                        'Visible',   'on');                                       
                    guidata(MLData.feature_display_handles{l},handles);
                else
                    MLData.feature_display_handles{l} = [];                    
                    MLData.feature_display_handles{l} = ...
                        plot(MLData.feature_display_axes{l}, MLData.Features{l}.feature(MLData.idx:MLData.idx+MLData.windowSize,:))';
                    for k = 1:length(MLData.feature_display_handles{l}),
                        set(MLData.feature_display_handles{l}(k),'Parent',MLData.feature_display_axes{l})                           
                        guidata(MLData.feature_display_handles{l}(k),handles);
                    end
                    % Very Quick Fix - Don't leave it like this
                    MLData.feature_display_centerlines{l} = line([1,1], ...
                        ylim(MLData.feature_display_axes{l}),...
                        'Color','r',...
                        'LineStyle','--',...
                        'LineWidth',1,...
                        'Parent', MLData.feature_display_axes{l});
                    guidata(MLData.feature_display_centerlines{l},handles);
                end
                
                set(MLData.feature_display_axes{l},'YDir','normal')                           
            
            case 'imagesc'
                MLData.feature_display_handles{l} = ...
                    imagesc(MLData.Features{l}.feature(MLData.idx:MLData.idx+MLData.windowSize,:)',...
                    'Parent',MLData.feature_display_axes{l});
                MLData.feature_display_centerlines{l} = line([1,1], ...
                    ylim(MLData.feature_display_axes{l}),...
                    'Color','r',...
                    'LineStyle','--',...
                    'LineWidth',1,...
                    'Parent', MLData.feature_display_axes{l});
                guidata(MLData.feature_display_handles{l},handles);
        end
    end
    set(MLData.feature_display_panels{l},'Position',display_panel_position(ssl,:))
    ssl = ssl+1;
end
unselected_features = setdiff(1:length(MLData.feature_display_axes),selected_features);
set(cell2mat(MLData.feature_display_panels(unselected_features)),'Visible','off');
set(cell2mat(MLData.feature_display_panels(selected_features)),'Visible','on');
set(cell2mat(MLData.feature_display_slider(unselected_features)),'Visible','off');
set(cell2mat(MLData.feature_display_slider(selected_features)),'Visible','on');
set(cell2mat(MLData.feature_display_axes(unselected_features)),'Visible','off');
set(cell2mat(MLData.feature_display_axes(selected_features)),'Visible','on');
set(cell2mat(MLData.feature_display_handles(selected_features)),'Visible','on');
set(cell2mat(MLData.feature_display_handles(unselected_features)),'Visible','off');
set(cell2mat(MLData.feature_display_centerlines(unselected_features)),'Visible','off');
set(cell2mat(MLData.feature_display_centerlines(selected_features)),'Visible','on');
set(cell2mat(MLData.feature_display_YLimCB(unselected_features)),'Visible','off');
set(cell2mat(MLData.feature_display_YLimCB(selected_features)),'Visible','on');

setappdata(handles.MTABrowser,'MLData',MLData);
guidata(hObject, handles);
updateCircle(hObject,handles,0)

function MLauxdataUpdatePlot(hObject,eventdata,handles)
MLData = getappdata(handles.MTABrowser,'MLData');
selectedFeatureName = get(handles.MLFnewFeatureName,'String');
selectedFeatureIndex =find(cellfun(@isequal,MLData.featurelabels,repmat({selectedFeatureName},1,MLData.num_of_features)),1);

if ~isempty(selectedFeatureIndex)
    try,
        delete(MLData.feature_display_handles{selectedFeatureIndex}),
        delete(MLData.feature_display_centerlines{selectedFeatureIndex}),
    end
    MLData.feature_display_handles{selectedFeatureIndex} = {0};
    MLData.feature_display_centerlines{selectedFeatureIndex} = {0};
    
    switch MLData.Features{selectedFeatureIndex}.plot_type,
        case 'line'            
            if size(MLData.Features{selectedFeatureIndex}.feature,2)==1,
                MLData.feature_display_handles{selectedFeatureIndex} =line('Parent', MLData.feature_display_axes{selectedFeatureIndex},...
                    'XData',     1:MLData.windowSize+1,...
                    'YData',     MLData.Features{selectedFeatureIndex}.feature(MLData.idx:MLData.idx+MLData.windowSize),...
                    'LineWidth', MLData.state_line_options.width, ...
                    'EraseMode', MLData.state_line_options.erase,...
                    'LineStyle', MLData.state_line_options.line_style,...
                    'Visible',   'on');
                guidata(MLData.feature_display_handles{selectedFeatureIndex},handles);                
            else                
                    MLData.feature_display_handles{selectedFeatureIndex} = [];                    
                    MLData.feature_display_handles{selectedFeatureIndex} = ...
                        plot(MLData.feature_display_axes{selectedFeatureIndex}, MLData.Features{selectedFeatureIndex}.feature(MLData.idx:MLData.idx+MLData.windowSize,:))';
                    for k = 1:length(MLData.feature_display_handles{selectedFeatureIndex}),
                        set(MLData.feature_display_handles{selectedFeatureIndex}(k),'Parent',MLData.feature_display_axes{selectedFeatureIndex})                           
                        guidata(MLData.feature_display_handles{selectedFeatureIndex}(k),handles);
                    end
                    % Very Quick Fix - Don't leave it like this
                    MLData.feature_display_centerlines{selectedFeatureIndex} = line([1,1], ...
                        ylim(MLData.feature_display_axes{selectedFeatureIndex}),...
                        'Color','r',...
                        'LineStyle','--',...
                        'LineWidth',1,...
                        'Parent', MLData.feature_display_axes{selectedFeatureIndex});
                    guidata(MLData.feature_display_centerlines{selectedFeatureIndex},handles);
              
            end
            set(MLData.feature_display_axes{selectedFeatureIndex},'YDir','normal')
            
        case 'imagesc'
            MLData.feature_display_handles{selectedFeatureIndex} = ...
                imagesc(MLData.Features{selectedFeatureIndex}.feature(MLData.idx:MLData.idx+MLData.windowSize,:)',...
                'Parent',MLData.feature_display_axes{selectedFeatureIndex});
            guidata(MLData.feature_display_handles{selectedFeatureIndex},handles);
            set(MLData.feature_display_axes{selectedFeatureIndex},'YDir','normal')
    end
    MLData.feature_display_centerlines{selectedFeatureIndex} = line([1,1], ...
        ylim(MLData.feature_display_axes{selectedFeatureIndex}),...
        'Color','r',...
        'LineStyle','--',...
        'LineWidth',1,...
        'Parent', MLData.feature_display_axes{selectedFeatureIndex});
    guidata(MLData.feature_display_centerlines{selectedFeatureIndex},handles);
end

setappdata(handles.MTABrowser,'MLData',MLData);
guidata(hObject,handles)
        
function MLauxdataDisplayOption_Callback(hObject, eventdata, handles)
MLData = getappdata(handles.MTABrowser,'MLData');
Tag = get(hObject,'Tag');
Style = get(hObject,'Style');
switch Style,
    
    case 'slider',
        switch Tag,
            case 'feature_display_slider',
                
                Value = get(hObject,'Value');
                index = find(cellfun(@isequal,MLData.feature_display_slider,repmat({hObject},1,length(MLData.feature_display_slider))),1);
                switch MLData.Features{index}.plot_type,
                    case 'line'
                        set(MLData.feature_display_axes{index},'YLim',[mean(MLData.Features{index}.feature(:))-Value,mean(MLData.Features{index}.feature(:))+Value]);
                    case 'imagesc'
                        caxis(MLData.feature_display_axes{index},[mean(MLData.Features{index}.feature(:))-Value,mean(MLData.Features{index}.feature(:))+Value]);
                end
                guidata(MLData.feature_display_axes{index},handles);
                setappdata(handles.MTABrowser,'MLData',MLData);
                guidata(hObject,handles);
        end
        
    case 'checkbox',
        switch Tag,         
            
            case 'feature_display_YLimCB',
                Value = get(hObject,'Value');
                index = find(cellfun(@isequal,MLData.feature_display_YLimCB,repmat({hObject},1,length(MLData.feature_display_YLimCB))),1);     
                if ~Value,
                    %ylim(MLData.feature_display_axes{index},'auto');
                    set(MLData.feature_display_axes{index},'YLimMode','auto');
                    set(MLData.feature_display_slider{index},'Visible','off');
                    %addlistener(MLData.feature_display_handles{index}, 'YData', 'PostSet', @(src, evnt) set(evnt.AffectedObject.Parent, 'YLim', [min(evnt.AffectedObject.YData) max(evnt.AffectedObject.YData)]));
                    guidata(MLData.feature_display_axes{index},handles);
                else
                    %ylim(MLData.feature_display_axes{index},'manual');
                    set(MLData.feature_display_axes{index},'YLimMode','manual');
                    set(MLData.feature_display_slider{index},'Visible','on');
                    guidata(MLData.feature_display_axes{index},handles);
                end
        end
end
setappdata(handles.MTABrowser,'MLData',MLData);
guidata(hObject,handles);

function MLFfeatures_loadTable(hObject, eventdata, handles)
MLData = getappdata(handles.MTABrowser,'MLData');
columnname =   {'Feature','visible','save feature'};
columnformat = {'char'   ,'logical','logical'};
columneditable=[false    , true    ,true];
featuresTableData = {};
for i = 1:length(MLData.Features),
    featuresTableData{i,1} = MLData.featurelabels{i};
    featuresTableData{i,2} = logical(MLData.selected_features(i));
    featuresTableData{i,3} = logical(MLData.Features{i}.ifsave);
end
set(handles.MLFfeatureTable,'ColumnName',columnname);
set(handles.MLFfeatureTable,'ColumnFormat',columnformat);
set(handles.MLFfeatureTable,'ColumnEditable',columneditable);
set(handles.MLFfeatureTable,'RowName',[]);
set(handles.MLFfeatureTable,'Data',featuresTableData);
set(handles.MLFfeatureTable,'ColumnWidth',{147,80,80,150})

set(handles.MLFnewFeatureName,'String',featuresTableData{1,1});
set(handles.MLFnewFeatureExpression,'String',MLData.Features{1}.expression);

guidata(handles.MLFfeatureTable, handles);
function MLFfeatureTable_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to MLFfeatureTable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
MLData = getappdata(handles.MTABrowser,'MLData');
col = eventdata.Indices(2);
row = eventdata.Indices(1);
switch col
    case 2
        MLData.selected_features(row) = ~MLData.selected_features(row);
    case 3
        MLData.Features{row}.ifsave = ~MLData.Features{row}.ifsave;
    otherwise
        return
end
setappdata(handles.MTABrowser,'MLData',MLData);
guidata(hObject, handles);
MLfeatureView_refresh(hObject,eventdata,handles);
function MLFfeatureTable_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to MLFfeatureTable (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
if ~isempty(eventdata.Indices);
    MLData = getappdata(handles.MTABrowser,'MLData');
    row = eventdata.Indices(1);
    featureTabelData = get(handles.MLFfeatureTable,'Data');
    set(handles.MLFnewFeatureName,'String',featureTabelData{row,1});
    set(handles.MLFnewFeatureExpression,'String',MLData.Features{row}.expression);
    displayTypes = cellstr(get(handles.MLFdisplayType,'String'));
    displayTypeIndex =find(cellfun(@isequal,displayTypes',repmat({MLData.Features{row}.plot_type},1,length(displayTypes))),1);
    set(handles.MLFdisplayType,'Value',displayTypeIndex);
    guidata(hObject, handles);
end


function MLFaddFeature_Callback(hObject, eventdata, handles)
% Get the name and expression of the feature which will be added
MLData = getappdata(handles.MTABrowser,'MLData');
newFeatureName = get(handles.MLFnewFeatureName,'String');
if isempty(newFeatureName)|sum(~cellfun(@isempty,regexpi(MLData.featurelabels,newFeatureName))),
    warndlg({'New feature name field is empty','or the name is arleady set','for another feature'}),
    return
end
newFeatureExpression = get(handles.MLFnewFeatureExpression,'String');
if isempty(newFeatureExpression),
    warndlg({'New feature expression field is empty'}),
    return
end

displayTypeIndex = get(handles.MLFdisplayType,'Value');
if displayTypeIndex == 1,
    warndlg({'Please select the feature''s display type'}),
    return
end

displayTypes = cellstr(get(handles.MLFdisplayType,'String'));
newFeaturePlotType = displayTypes{displayTypeIndex};

MLData.featurelabels{end+1} = newFeatureName;

MLData.Features{end+1} = MTAFeature(newFeatureName,newFeatureExpression,0,newFeaturePlotType);
MLData.num_of_features = MLData.num_of_features + 1;
MLData.selected_features(end+1) = false;
MLData.feature_display_handles(end+1)     ={0};
MLData.feature_display_YLimCB(end+1)      ={0};
MLData.feature_display_centerlines(end+1) ={0};
MLData.feature_display_axes(end+1)        ={0};
MLData.feature_display_slider(end+1)      ={0};
MLData.feature_display_panels(end+1)      ={0};

setappdata(handles.MTABrowser,'MLData',MLData);
guidata(hObject, handles);
MLFfeatures_loadTable(hObject, eventdata, handles);
MLfeatureView_refresh(hObject,eventdata,handles);
function MLFremoveFeature_Callback(hObject, eventdata, handles)
% Get the name and expression of the feature which will be added
MLData = getappdata(handles.MTABrowser,'MLData');
newFeatureName = get(handles.MLFnewFeatureName,'String');
removeFeatureIndex = find(~cellfun(@isempty,regexpi(MLData.featurelabels,newFeatureName)));
keptFeatures = setdiff(1:MLData.num_of_features,removeFeatureIndex);
MLData.featurelabels = MLData.featurelabels(keptFeatures);
MLData.Features = MLData.Features(keptFeatures);
MLData.num_of_features = MLData.num_of_features - 1;
MLData.selected_features = MLData.selected_features(keptFeatures);


if ~isequal(MLData.feature_display_handles{removeFeatureIndex},0)|...
    ishandle(MLData.feature_display_handles{removeFeatureIndex});
    delete(MLData.feature_display_handles{removeFeatureIndex});
    delete(MLData.feature_display_centerlines{removeFeatureIndex});
    delete(MLData.feature_display_axes{removeFeatureIndex});
    delete(MLData.feature_display_slider{removeFeatureIndex});
    delete(MLData.feature_display_YLimCB{removeFeatureIndex});
    delete(MLData.feature_display_panels{removeFeatureIndex});
end

if isempty(keptFeatures),
    MLData.feature_display_handles     ={};
    MLData.feature_display_centerlines ={};
    MLData.feature_display_axes        ={};
    MLData.feature_display_panels      ={};
    MLData.feature_display_YLimCB      ={};
    MLData.feature_display_slider      ={};
else
    MLData.feature_display_handles     = MLData.feature_display_handles(keptFeatures);
    MLData.feature_display_YLimCB     = MLData.feature_display_YLimCB(keptFeatures);
    MLData.feature_display_centerlines = MLData.feature_display_centerlines(keptFeatures);
    MLData.feature_display_axes        = MLData.feature_display_axes(keptFeatures);
    MLData.feature_display_slider        = MLData.feature_display_slider(keptFeatures);
    MLData.feature_display_panels      = MLData.feature_display_panels(keptFeatures);
end
setappdata(handles.MTABrowser,'MLData',MLData);
guidata(hObject, handles);
MLFfeatures_loadTable(hObject, eventdata, handles);
MLfeatureView_refresh(hObject,eventdata,handles);
function MLFmoveFeatureUp_Callback(hObject, eventdata, handles)
function MLFmoveFeatureDown_Callback(hObject, eventdata, handles)
function MLFnewFeatureName_Callback(hObject, eventdata, handles)
function MLFdisplayType_Callback(hObject, eventdata, handles)
MLData = getappdata(handles.MTABrowser,'MLData');
displayTypeIndex = get(hObject,'Value');
if displayTypeIndex == 1,
    warndlg({'Please select the feature''s display type'}),
    return
end
displayTypes = cellstr(get(hObject,'String'));
newFeaturePlotType = displayTypes{displayTypeIndex};
selectedFeatureName = get(handles.MLFnewFeatureName,'String');
selectedFeatureIndex =find(cellfun(@isequal,MLData.featurelabels,repmat({selectedFeatureName},1,MLData.num_of_features)),1);
if isempty(selectedFeatureIndex)|displayTypeIndex==1,
    return
else    
    MLData.Features{selectedFeatureIndex}.plot_type = newFeaturePlotType;
end
setappdata(handles.MTABrowser,'MLData',MLData);
MLauxdataUpdatePlot(hObject,eventdata,handles)
guidata(hObject,handles);


function MLFnewFeatureExpression_Callback(hObject, eventdata, handles)
function MLFfetOpen_Callback(hObject, eventdata, handles)
Session = getappdata(handles.MTABrowser,'Session');
MLData = getappdata(handles.MTABrowser,'MLData');
if isempty(Session.fet),
    fet_name = 'default';
else
    fet_name = Session.fet.mode;
end

[filename,filepath] = uigetfile(['*.' Session.trialName '.fet.*'],'Open Feature Collection:',...
                     fullfile(Session.spath, [Session.filebase '.fet.' fet_name '.mat']),...
                     'MultiSelect','off');
                 
if ~ischar(filename)||~ischar(filepath),                 
    return
else
    load([filepath filename]);
    Session.fet = Fet;
    Features = Session.fet.Features;
    if strcmp(Session.fet.mode,'default'),        
        Session.fet = MTAFet(Session,fet_name,1);
        featurelabels = {'default'};
        Features{1} = MTAFeature('default','zeros(Session.xyz.size(1),1);',0,'line');
        Session.fet.Features = MLData.Features;
    else
        if ~isempty(Features)
            featurelabels = cell(1,length(Features));
            for i = 1:length(Features)
                featurelabels{i} = Features{i}.label;
            end
        end
    end
    MLData.Features = Features;
    MLData.featurelabels = featurelabels;
    MLData.num_of_features = length(Features);
    MLData.selected_features = false(1,MLData.num_of_features);

    setappdata(handles.MTABrowser,'MLData',MLData)
    set(handles.MLFtitle,'String',['Bhavioral Features - ' filename]);
    guidata(hObject, handles);
    MLFfeatures_loadTable(hObject, eventdata, handles);
    MLfeatureView_refresh(hObject,eventdata,handles);
end
function MLFfetSave_Callback(hObject, eventdata, handles)
Session = getappdata(handles.MTABrowser,'Session');
MLData = getappdata(handles.MTABrowser,'MLData');
if isempty(Session.fet),
    fet_name = 'default';
else
    fet_name = Session.fet.mode;
end
[filename,filepath] = uiputfile(['*' Session.trialName '.fet.*'],'Open Feature Collection:',...
                     fullfile(Session.spath, [Session.filebase '.fet.' fet_name '.mat']));

if ~ischar(filename)||~ischar(filepath), 
    return
else
    points = regexp(filename,'[.]');
    fileparts = [1 points+1; points-1 length(filename)]';
    if size(fileparts,1)==6,
        fet_name = filename(fileparts(5,1):fileparts(5,2));
    else
        error('Fet filename does not match standard format,')
    end
    if ~exist([filepath filename],'file'),
        Fet = MTAFet(Session,fet_name,1);
    else
        load([filepath filename]);
    end
    Fet.Features = MLData.Features;
    Fet.save(Session,1);
    Session.fet = Fet;
    Session.save();
end

function MLFdisplayType_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MLFdisplayType (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function MLFnewFeatureName_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function MLFnewFeatureExpression_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function MLFaddFeature_CreateFcn(hObject, eventdata, handles)
icon = imread('add.png');
set(hObject,'Cdata',icon);
guidata(hObject, handles);
function MLFremoveFeature_CreateFcn(hObject, eventdata, handles)
icon = imread('subtract.png');
set(hObject,'Cdata',icon);
guidata(hObject, handles);
function MLFmoveFeatureUp_CreateFcn(hObject, eventdata, handles)
icon = imread('arrow_up.png');
set(hObject,'Cdata',icon);
guidata(hObject, handles);
function MLFmoveFeatureDown_CreateFcn(hObject, eventdata, handles)
icon = imread('arrow_down.png');
set(hObject,'CData',icon);
guidata(hObject, handles);
function MLFfetOpen_CreateFcn(hObject, eventdata, handles)
icon = imread('open.png');
set(hObject,'Cdata',icon);
guidata(hObject, handles);
function MLFfetSave_CreateFcn(hObject, eventdata, handles)
icon = imread('save.png');
set(hObject,'Cdata',icon);
guidata(hObject, handles);
%% MLF feature Functions end


function MLviewoptions_Callback(hObject, eventdata, handles)
AppStateChange(hObject,eventdata,handles);

function MLkeyState(hObject,eventdata,handles)
%% MLKEYSTATE
whatkey = double(get(handles.MTABrowser,'CurrentCharacter'));
if ~whatkey,return,end
switch whatkey
    case double(' ')
        setappdata(handles.MLapp,'paused',~getappdata(handles.MLapp,'paused'));
    case double('e')
        MLstateErase_Callback(hObject, eventdata, handles)
    case 29
        MLplayforward_Callback(hObject, eventdata, handles)
    case 28
        MLplayreverse_Callback(hObject, eventdata, handles)
    otherwise
        MLData = getappdata(handles.MTABrowser,'MLData');
        if MLData.current_label == strfind(MLData.keys(MLData.selected_states),char(whatkey));
            MLData.flag_tag=~MLData.flag_tag;
        else
            if isempty(MLData.current_label),
                MLData.current_label = 0;
            else
                MLData.current_label = strfind(MLData.keys(MLData.selected_states),char(whatkey));
                MLData.flag_tag=1;
            end
        end
        if MLData.flag_tag
            set(handles.MLstateLabel,'BackgroundColor',[0 1 1])
            MLData.erase_tag=false;
            set(handles.MLstateErase,'BackgroundColor',[0.941176 0.941176 0.941176])
        else
            set(handles.MLstateLabel,'BackgroundColor',[0.941176 0.941176 0.941176])
        end
        setappdata(handles.MTABrowser,'MLData',MLData);
end
set(handles.MTABrowser,'CurrentCharacter','0');
guidata(handles.MTABrowser,handles);
function updateCircle(hObject,handles,idx)
%% UpdateCircle

play_speed = getappdata(handles.MLapp,'play_speed');
pause(0.05/play_speed);
if ~getappdata(handles.MLapp,'paused'),
    MLData = getappdata(handles.MTABrowser,'MLData');

    skip = round(play_speed/10);
    
    selected_states = find(MLData.selected_states);
    % Label state
    if MLData.flag_tag&MLData.current_label,
        MLData.States{selected_states(MLData.current_label)}.data(MLData.idx:MLData.idx+skip*idx+idx) = MLData.current_label/2;
    end
    % Erase state
    if MLData.erase_tag&MLData.current_label,
        MLData.States{selected_states(MLData.current_label)}.data(MLData.idx:MLData.idx+skip*idx+idx) = 0;
    end

    MLData.idx = MLData.idx+skip*idx+idx;

    if ( MLData.idx < 1 ), 
        MLData.idx = MLData.sessionLength-MLData.windowSize; 
    elseif ( MLData.idx > MLData.sessionLength-MLData.windowSize ),
        MLData.idx = 1;
    end

    set(handles.MLposition_slider,'Value', MLData.idx)
    for kk=1:MLData.num_of_sticks,
        set(MLData.sticks{kk},'XData',[MLData.xyzpos(MLData.idx,MLData.markerConnections(kk,1),1),MLData.xyzpos(MLData.idx,MLData.markerConnections(kk,2),1)],...
                               'YData',[MLData.xyzpos(MLData.idx,MLData.markerConnections(kk,1),2),MLData.xyzpos(MLData.idx,MLData.markerConnections(kk,2),2)],...
                               'ZData',[MLData.xyzpos(MLData.idx,MLData.markerConnections(kk,1),3),MLData.xyzpos(MLData.idx,MLData.markerConnections(kk,2),3)])
    end
    for kk=1:MLData.num_of_markers,
        set(MLData.markers{kk},'XData',[MLData.xyzpos(MLData.idx,kk,1),MLData.xyzpos(MLData.idx,kk,1)],...
                               'YData',[MLData.xyzpos(MLData.idx,kk,2),MLData.xyzpos(MLData.idx,kk,2)],...
                               'ZData',[MLData.xyzpos(MLData.idx,kk,3),MLData.xyzpos(MLData.idx,kk,3)])
    end
    for kk=find(MLData.selected_states),
        set(MLData.state_lines{kk},'XData',MLData.index_range(MLData.idx,1):MLData.index_range(MLData.idx,2),...
                                   'YData',MLData.States{kk}.data(MLData.index_range(MLData.idx,1):MLData.index_range(MLData.idx,2)),...
                                   'LineStyle','-')
    end
    for kk=find(MLData.selected_features),
        switch MLData.Features{kk}.plot_type,
            case 'line'
                if length(MLData.feature_display_handles{kk})==1,                
                set(MLData.feature_display_handles{kk},...
                    'XData',MLData.index_range(MLData.idx,1):MLData.index_range(MLData.idx,2),...
                    'YData',MLData.Features{kk}.feature(MLData.index_range(MLData.idx,1):MLData.index_range(MLData.idx,2)),...
                    'LineStyle','-')
                else
                    for l = 1:length(MLData.feature_display_handles{kk}),
                        set(MLData.feature_display_handles{kk}(l),...
                            'XData',MLData.index_range(MLData.idx,1):MLData.index_range(MLData.idx,2),...
                            'YData',MLData.Features{kk}.feature(MLData.index_range(MLData.idx,1):MLData.index_range(MLData.idx,2),l),...
                            'LineStyle','-')
                    end
                end
                
                set(MLData.feature_display_centerlines{kk},'XData',[MLData.idx,MLData.idx],...
                    'YData',[min(min(MLData.Features{kk}.feature(MLData.index_range(MLData.idx,1):MLData.index_range(MLData.idx,2))))-1,max(max(MLData.Features{kk}.feature(MLData.index_range(MLData.idx,1):MLData.index_range(MLData.idx,2))))+1]);
                
            case 'imagesc'
                set(MLData.feature_display_handles{kk},...
                    'XData',[MLData.index_range(MLData.idx,1):MLData.index_range(MLData.idx,2)]',...
                    'YData',1:size(MLData.Features{kk}.feature,2),...
                    'CData',MLData.Features{kk}.feature(MLData.index_range(MLData.idx,1):MLData.index_range(MLData.idx,2),:)');
                set(MLData.feature_display_centerlines{kk},'XData',[MLData.idx,MLData.idx],...
                    'YData',[1,size(MLData.Features{kk}.feature,2)])

        end
    end
    
    set(cell2mat(MLData.feature_display_axes(MLData.selected_features)),'XLim',[MLData.idx-MLData.windowSize, MLData.idx+MLData.windowSize])
    set(handles.MLstateView,'XLim',[MLData.idx-MLData.windowSize, MLData.idx+MLData.windowSize],'YLim',[MLData.statesRange(1) MLData.statesRange(2)])
    set(MLData.state_centerline,'Xdata',[MLData.idx,MLData.idx] )
    setappdata(handles.MTABrowser,'MLData',MLData);
    set(handles.MLindexinput,'String',num2str(MLData.idx));
    set(handles.MLtimeinput,'String',num2str(MLData.t(MLData.idx)));
    drawnow
    if MLData.recording,
        MLData.record{end+1} = getframe(handles.MLxyzView);
        MLData.record{end}.index = MLData.idx;        
        setappdata(handles.MTABrowser,'MLData',MLData);
    end
    
end
function MLstateLabel_Callback(hObject, eventdata, handles)
% hObject    handle to MLstateLabel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

MLData = getappdata(handles.MTABrowser,'MLData');
MLData.flag_tag = ~MLData.flag_tag;

if MLData.flag_tag
    set(handles.MLstateLabel,'BackgroundColor',[0 1 1])
    MLData.erase_tag=false;
    set(handles.MLstateErase,'BackgroundColor',[0.941176 0.941176 0.941176])
else
    set(handles.MLstateLabel,'BackgroundColor',[0.941176 0.941176 0.941176])
end

setappdata(handles.MTABrowser,'MLData',MLData);
guidata(hObject,handles);
function MLstateErase_Callback(hObject, eventdata, handles)
% hObject    handle to MLstateErase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
MLData = getappdata(handles.MTABrowser,'MLData');
MLData.erase_tag = ~MLData.erase_tag;

if MLData.erase_tag
    set(handles.MLstateErase,'BackgroundColor',[0 1 1])
    MLData.flag_tag=false;
    set(handles.MLstateLabel,'BackgroundColor',[0.941176 0.941176 0.941176])
else
    set(handles.MLstateErase,'BackgroundColor',[0.941176 0.941176 0.941176])
end

setappdata(handles.MTABrowser,'MLData',MLData);
guidata(hObject,handles);

%% MLapp End






%% LSapp Start
% --- Executes on button press in BSlfpStates.
function BSlfpStates_Callback(hObject, eventdata, handles)
% hObject    handle to BSlfpStates (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%% Set Browser State
BStateChange(hObject,eventdata,handles);
% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double
% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double
% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% 
% % --- Executes on button press in LSstates.
% function togglebutton8_Callback(hObject, eventdata, handles)
% % hObject    handle to LSstates (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: get(hObject,'Value') returns toggle state of LSstates
% 
% 
% % --- Executes on button press in LSfeatures.
% function togglebutton9_Callback(hObject, eventdata, handles)
% % hObject    handle to LSfeatures (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hint: get(hObject,'Value') returns toggle state of LSfeatures
% 

%       {struct('object',handles.MLSapp,...
%               'properties',struct(...
%                   'Position',[],...
%                   'Visible' ,'on')),...
%        struct('object',handles.MLFapp,...
%               'properties',struct(...
%                   'Position',[],...
%                   'Visible' ,'on')),...
%        struct('object',handles.MLVOapp,...
%               'properties',struct(...
%                   'Position',[],
%                   'Visible' ,'on')),...

%       {struct('object', handles.MLSapp,...
%               'properties',struct(...
%                   'Position',[0.1,0.1,0.1,0.1],...
%                   'Visible' ,'off')),...
%        struct('object', handles.MLFapp,...
%               'properties',struct(...
%                   'Position',[0.1,0.1,0.1,0.1],...
%                   'Visible' ,'off')),...
%        struct('object', handles.MLVOapp,...
%               'properties',struct(...
%                   'Position',[0.1,0.1,0.1,0.1],...
%                   'Visible' ,'off')),...


function LSveiwoptions_Callback(hObject, eventdata, handles)

function LSstates_Callback(hObject, eventdata, handles)

function LSfeatures_Callback(hObject, eventdata, handles)
