function [consistent,oriEvtTotalTime,adjEvtTotalTime] = check_xyz_nlx_consistency(Session)

if ~isa(Session,'MTASession'),
    Session = MTASession(Session);
end

%% Get total recording time for concatinated nlx recordings
adjEvt = LoadEvents([Session.spath.nlx Session.name '.all.evt']);


re = 'Starting Recording';
adjEvtStart = adjEvt.time(~cellfun(@isempty,regexp(adjEvt.description,re)));
re = 'Stopping Recording';
adjEvtStop = adjEvt.time(~cellfun(@isempty,regexp(adjEvt.description,re)));


if isempty(adjEvtStart), error('adjusted start event is missing'),end
if isempty(adjEvtStop),  error('adjusted stop event is missing'),end

if sum(size(adjEvtStart))~=sum(size(adjEvtStop)),
    error(['number of start and stop events do not match in adjusted evt file']);
end

adjEvtTotalTime = sum(adjEvtStop-adjEvtStart);

%% Find original event files and get total recording time for each
%% sub session
files = dir([Session.spath.nlx]);
re = ['nlx.evt'];
oriEvtFileList = {files(~cellfun(@isempty,regexp({files.name},re))).name};

oriEvtTotalTime = 0;
for file = 1:length(oriEvtFileList),
    oriEvt = LoadEvents([Session.spath.nlx oriEvtFileList{file}]);

    re = 'Starting Recording';
    oriEvtStart = oriEvt.time(~cellfun(@isempty,regexp(oriEvt.description,re)));
    re = 'Stopping Recording';
    oriEvtStop = oriEvt.time(~cellfun(@isempty,regexp(oriEvt.description,re)));


    if isempty(oriEvtStart), error('Original start event is missing'),end
    if isempty(oriEvtStop),  error('Original stop event is missing'),end


    if sum(size(oriEvtStart))~=sum(size(oriEvtStop)),
        error(['number of start and stop events do not match in ' ...
               'original evt file: ' Session.spath.nlx oriEvtFileList{file}])
    end

    oriEvtTotalTime = sum(oriEvtStop-oriEvtStart)+oriEvtTotalTime;
end

if round(oriEvtTotalTime)==round(adjEvtTotalTime),
    consistent = 1;
else
    consistent = 0;
end

