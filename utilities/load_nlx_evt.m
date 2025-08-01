function events = load_nlx_evt(filename, varargin)
% load_nlx_evt - Read events from file.
%
%  USAGE
%
%    events = LoadEvents(filename)
%
%    filename            event file name
%
%  OUTPUT 
% events = 
%
%  struct with fields:
%
%           time: [196x1 double]
%    description: {196x1 cell}
%         Labels: {10x1 cell}
%            Clu: [196x1 double]
%
%
% Copyright (C) 2004-2006 by Michaël Zugaro
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
%
% Repurposed by Justin Graboski


[LabelNumeric] = DefaultArgs(varargin,{0});

events.time = [];
events.description = [];
if ~exist(filename),
	error(['File ''' filename ''' not found.']);
end

file = fopen(filename,'r');
if file == -1,
	error(['Cannot read ' filename ' (insufficient access rights?).']);
end

while ~feof(file),
	time = fscanf(file,'%f',1);
	if isempty(time),
		if feof(file), break; end
		error(['Failed to read events from ' filename ' (possibly an empty file).']);
	end
	events.time(end+1,1) = time;
	line = fgetl(file);
	start = regexp(line,'[^\s]','once');
	events.description{end+1,1} = sscanf(line(start:end),'%c');
end
fclose(file);

% Convert to seconds
if ~isempty(events.time), events.time(:,1) = events.time(:,1) / 1000; end
% make numeric labels

[events.labels, ~, events.Clu] = unique(events.description);
if LabelNumeric
% check for the order
    CurOrder = cellfun(@str2num, events.labels);
    [~, NewOrder] = sort(CurOrder);
    events.labels = events.labels(NewOrder);
    events.clu = CurOrder(events.clu);    
end



