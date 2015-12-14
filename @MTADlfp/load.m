function Data = load(Data,Session,varargin)
%load(Data,Session,varargin)
%[channels,gselect,periods] = DefaultArgs(varargin,{[],{'AnatGrps',1,1},[]});
[channels,gselect,periods] = DefaultArgs(varargin,...
                                         {[],{'AnatGrps',1,1},[]});
if isempty(periods),
    if Session.sync.sampleRate~=1,Session.sync.resample(1);end
    periods = round(Session.sync([1,end]).*Session.lfp.sampleRate);
elseif isa(periods,'MTADepoch'),
    periods.resample(Session.lfp.sampleRate);
    periods = periods.data + Session.sync(1)*Session.lfp.sampleRate;
end
Par = LoadPar(fullfile(Session.spath, [Session.name '.xml']));
if isempty(channels)
    channels = Par.(gselect{1})(gselect{2}).Channels(gselect{3});
end
Data.data = LoadBinary(Data.fpath,channels,...
                       Par.nChannels,[],[],[],periods)';
Data.data(Data.data==0)=1;
%Session.resync(Data);
end
