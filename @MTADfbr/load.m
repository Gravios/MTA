function Data = load(Data,Session,varargin)
%load(Data,Session,varargin)
%[channels,gselect,periods] = DefaultArgs(varargin,{[],{'AnatGrps',1,1},[]});
[channels,gselect,periods] = DefaultArgs(varargin,...
                                         {[],{'AnatGrps',1,1},[]});
if isempty(periods),
    if Session.sync.sampleRate~=1,Session.sync.resample(1);end
    periods = round(Session.sync([1,end]).*Session.fbr.sampleRate);
elseif isa(periods,'MTADepoch'),
    periods.resample(Session.fbr.sampleRate);
    periods = periods.data + Session.sync(1)*Session.fbr.sampleRate;
end
Par = LoadPar(fullfile(Session.spath, [Session.name '.xml']));
if isempty(channels)
    channels = ':';
end
Data.clear();

% LOAD original data ... all of it ... I know it's a waste of memory
load(fullfile(Session.spath,[Session.name,'.fbr']),'-mat');
for period = periods'
    Data.data = cat(1,Data.data,fbrData(period(1):period(2),channels));
end


end
