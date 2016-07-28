function evts = convert_nlx_events_to_mta_epochs(Trial,eventLabels,varargin)
% Use LoadEvents to import entries from the [Trial.name].*.evt file 
% and encapsulate in a MTADepoch object
%
% varargin:
%    eventFileType
%    label
%    key
%    name
%
% first case:
%    eventLabels is a 1x1 cell.
%    label (default) = eventLabels{1};
%    key   (default) = 'x';
%    name  (default) = '';
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
eventLabels = {'teleport'};
label = 'shifted';
key = 'x';

%% MAIN 
events = LoadEvents(fullfile(Trial.spath, [Trial.name '.' eventFileType '.evt']));

eClu = find(~cellfun(@isempty,regexp(events.Labels,eventLabels)));

eTime = events.time(events.Clu==eClu);


