function Data = load(Data,Session,varargin)
%load(Data,Session,varargin)
% channels - numericArray: select,periods] = DefaultArgs(varargin,{[],{'AnatGrps',1,1},[]});
%
% parameter file : xml file (ndm toolbox format) containing recording session info
% 
% TODO : switch loading of parameter file to lfp.model 
%
    
% DEFARGS ------------------------------------------------------------------------------------------
[channels,gselect,periods] = DefaultArgs(varargin,{[],{'AnatGrps',1,1},[]});
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

% REMOVE AFTER UPDATE ---------------------------------------------
% LOAD parameter file                                             %
par = LoadPar(fullfile(Session.spath, [Session.name '.xml']));    %
if isempty(channels)                                              %
    channels = par.(gselect{1})(gselect{2}).Channels(gselect{3}); %
end                                                               %
                                                                  %
switch Data.ext                                                   %
  case 'lfp'                                                      %
    Data.sampleRate = par.lfpSampleRate;                          %
  case 'dat'                                                      %
    Data.sampleRate = par.SampleRate;                             %
end                                                               %
%------------------------------------------------------------------

if isempty(periods),
    if Session.sync.sampleRate~=1,Session.sync.resample(1);end
    periods = round(Session.sync([1,end]).*Data.sampleRate);
elseif isa(periods,'MTADepoch'),
    periods.resample(Data.sampleRate);
    periods = periods.data + Session.sync(1)*Data.sampleRate;
else
    periods = round((Session.sync(1)+periods)*Data.sampleRate);
end


try,
    Data.data = LoadBinary(Data.fpath,channels,par.nChannels,[],[],[],periods)';
catch
    Data.data = LoadBinary(Data.fpath,channels,par.nChannels,[],[],[],periods)';
end

Data.data(Data.data==0)=1;
%Session.resync(Data);
end

% END MAIN -----------------------------------------------------------------------------------------