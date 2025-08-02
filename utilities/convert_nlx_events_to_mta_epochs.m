function events = convert_nlx_events_to_mta_epochs(Trial,varargin)
% Use LoadEvents to import entries from the [Trial.name].*.evt file 
% and encapsulate in a MTADepoch object
%
% varargin:
%    eventLabels
%    eventFileType
%    label
%    key
%    name
%
% first case:
%    eventLabels is a 1x1 cell.
%    eventLabels (default) = ''
%    label (default) = eventLabels{1};
%    key   (default) = 'x';
% 
%    single tag represents condition switch between to condition

%% Default Arguments
defargs = struct('eventFileType',       'all',     ...
                 'label',               'alt',     ...
                 'key',                 '
[eventFileType,label,key,name] = DefaultArgs(varargin,

%% Diagnostic Vars
MTAstartup('vr_exp');
Trial = MTATrial.validate('Ed10-20140820.rov.all');
eventLabels = 'teleport';
eventFileType = 'all';
label = 'shifted';
key = 'x';

%% MAIN 
events = LoadEvents(fullfile(Trial.spath, [Trial.name '.' eventFileType '.evt']));

eClu = find(~cellfun(@isempty,regexp(events.Labels,eventLabels)));

eTime = events.time(events.Clu==eClu);
eTime = round((eTime - Trial.sync.data(1)).*Trial.xyz.sampleRate);

if mod(numel(eTime),2)==0,
    eTime = [eTime,Trial.sync(end)];
end

events = MTADepoch(Trial.spath,...
                   Trial.filebase,...
                   reshape(eTime,2,[])',...
                   Trial.xyz.sampleRate,...
                   Trial.sync.copy,...
                   Trial.sync.data(1),...
                   [],[],[],label,key);
