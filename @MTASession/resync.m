function Data = resync(Session,Data,varargin)
% Data = resync(Session,Data,varargin)
% Resync uses the sync object of a session or one given through the
% varargin option. 
%
%  Session: MTASession, Data object holding all session information
%
%  Data:    MTAData,    Data object targeted for resynchronization.
%
%  sync:    MTADepoch,  A set of periods defining what data needs to
%                       be loaded. 
%           double,     Set of periods which have to be specifed in
%                       indecies in the sampling rate of the Data
%                       object
%
[sync] = DefaultArgs(varargin,{[]});

%global diagnostic;

%diagnostic,
%Data = Data.copy;

switch class(Data)
  case 'MTADepoch'
    if ~strcmp(Data.path,Session.spath),
        Data.updatePath(Session.spath);
    end
    
    if exist(Data.fpath,'file')&&~isempty(Data.filename),
        Data.load;
    end
end

if ~isempty(sync)
    switch class(sync)
      case 'double',
        msg.message = 'The provided synchronization periods are empty.';
        msg.identifier = 'MTASession:resync:EmptySync';
        msg.stack = dbstack;
        if isempty(sync),error(msg),end
        Data.sync.data = sync./Data.sampleRate+Session.sync.origin;
        
      case 'MTADepoch'
        Data.sync.sync = sync.copy;
    end
else
    Data.sync.sync = Session.sync.copy;
end




if isa(Data,'MTADepoch'),
    %%%%%%%%%%%%%%%%%%% REDO %%%%%%%%%%%%%%%%%%%
    %Data.resample(1);
% $$$        Data.data = IntersectRanges(Data.data+Data.origin,Data.sync.sync.data+Data.sync.sync.origin-1)-Data.sync.sync(1);
% $$$             Data.origin = Data.sync.sync(1)+1;
    if Data.origin ~= Data.sync.sync.data(1),
        indShift = round((Data.origin - ...
                          Data.sync.sync.data(1))*Data.sampleRate);            
    else
        indShift = 0;
    end
    syncp = Data.sync.sync.copy;
    syncp.resample(Data.sampleRate);
    
    Data.data = IntersectRanges(Data.data+indShift,syncp.data-syncp.data(1)+1);
    Data.origin = Data.sync.sync.data(1);

    
    return
    %%%%%%%%%%%%%%%%%%% REDO %%%%%%%%%%%%%%%%%%%

elseif isa(Data,'MTAData'),
    % TODO: Need to deal with lfp resync eventually.
    if isa(Data,'MTADlfp')||Data.isempty,
        Data.sync.sync = Session.sync.copy;
        return
    end
    
    % The periods when the data was recorded
    dataEpoch = Data.sync.copy;
    dataEpoch.cast('TimeSeries',Data.sampleRate,'absolute');
    dataOrigin = round(Data.sync(1)*Data.sampleRate);
    
    % The periods of data which are already loaded
    loadedData = ones(Data.size(1),1);
    try
        loadedData(Data.data(:,1,1,1,1)==0) = 0;
    end
    loadedData = cat(1,zeros(dataOrigin,1),loadedData);
    tailbuff = dataEpoch.size(1)-size(loadedData,1);
    if tailbuff<0,
        loadedData = cat(1,loadedData,zeros(dataEpoch.size(1)-size(loadedData,1),1));
    else
        loadedData = loadedData(1:end+tailbuff);
    end
    
    loadedDataEnd = find(loadedData==1,1,'last');
    
    % The desired synchronization periods
    syncEpoch = Data.sync.sync.copy;
    syncEpoch.cast('TimeSeries',Data.sampleRate,'absolute');
    try
        syncEpoch.data = syncEpoch.data(1:dataEpoch.size(1));
    catch err
        warning('Resync may be shifted 1 index');
        syncEpoch.data = [0;syncEpoch.data];
    end

    newOrigin = find(syncEpoch.data==1,1,'first');


    %Diagnostic
    %diagnostic,
% $$$                   figure,
% $$$                   plot(dataEpoch.data)
% $$$                   ylim([-3,3])
% $$$                   hold on
% $$$                   plot(loadedData-dataEpoch.data,'r')
% $$$                   plot(loadedData-syncEpoch.data,'g')
% $$$                   plot(dataEpoch.data-syncEpoch.data,'c')
% $$$                   plot(loadedData-syncEpoch.data+dataEpoch.data,'m')
% $$$ 
% $$$                   figure,
% $$$                   plot(dataEpoch.data)
% $$$                   ylim([-3,3])
% $$$                   hold on
% $$$                   plot(loadedData+.1,'r')
% $$$                   plot(syncEpoch.data+.2,'g')
% $$$                   plot(dataEpoch.data,'c')
% $$$                   plot(loadedData-syncEpoch.data+dataEpoch.data,'m')
% $$$ 
% $$$                   ds = load(Data.fpath)
% $$$                   figure,
% $$$                   plot(ds.data(:,1,1)==eps)
% $$$                   hold on
% $$$                   ,plot(ds.data(:,1,1)==0,'r')+.5
% $$$                   ylim([-3,3])
% $$$ 
% $$$                   figure
% $$$                   plot(Data.data(:,1,1)==eps)
% $$$                   hold on
% $$$                   plot((Data.data(:,1,1)==0)+.5,'r')
% $$$                   ylim([-3,3])
% $$$                   
    %%Trim ends
    endSync = syncEpoch.size(1);
    endShiftIndex = endSync - loadedDataEnd;
    endShiftIndex(endShiftIndex==0)=1;           
    startShiftIndex = newOrigin-dataOrigin;
    startShiftIndex(startShiftIndex==0)=1;           
    if endShiftIndex < 0,
        if startShiftIndex < 0,
            Data.data = cat(1,zeros([abs(startShiftIndex),Data.size(2:end)]),Data.data(1:abs(endSync-dataOrigin+1),:,:,:,:));
            dataEpoch.data = dataEpoch.data(newOrigin:endSync);
            loadedData = loadedData(newOrigin:endSync);
            syncEpoch.data = syncEpoch.data(newOrigin:endSync);
        else
            Data.data = Data.data(startShiftIndex:endSync-newOrigin,:,:,:,:);
            dataEpoch.data = dataEpoch.data(newOrigin:endSync);
            loadedData = loadedData(newOrigin:endSync);
            syncEpoch.data = syncEpoch.data(newOrigin:endSync);
        end
    else
        if endShiftIndex == 1, endShiftIndex = 0; end
        if startShiftIndex < 0,
            Data.data = cat(1,zeros([abs(startShiftIndex),Data.size(2:end)]),Data.data,zeros([abs(endShiftIndex),Data.size(2:end)]));
            dataEpoch.data = dataEpoch.data(newOrigin:endSync);
            loadedData = loadedData(newOrigin:endSync);
            syncEpoch.data = syncEpoch.data(newOrigin:endSync);
        else
            Data.data = cat(1,Data.data(startShiftIndex:end,:,:,:,:),zeros([abs(endShiftIndex),Data.size(2:end)]));
            dataEpoch.data = dataEpoch.data(newOrigin:endSync);
            loadedData = loadedData(newOrigin:endSync);
            syncEpoch.data = syncEpoch.data(newOrigin:endSync);
        end
    end
    
    %Diagnostic
    %diagnostic,
% $$$                  figure,
% $$$                  plot(loadedData)
% $$$                  ylim([-2,2])
% $$$                  hold on
% $$$                  plot(loadedData-dataEpoch.data,'r')
% $$$                  plot(loadedData-syncEpoch.data,'g')
% $$$                  plot(dataEpoch.data-syncEpoch.data,'c')
    %%keyboard
    
    syncFeature = (loadedData-dataEpoch.data).*syncEpoch.data-syncEpoch.data;
    syncDataIndex = syncFeature==-2;
    syncDataPeriods = ThreshCross(syncDataIndex,0.5,3);
    syncZeroIndex = syncFeature==0;
    
    if ~isempty(syncDataPeriods),
        keyboard
        %syncshift = 0;
        %syncshift = round(Data.sync(1).*Data.sampleRate)-newOrigin-1;
        %syncshift = Data.sync(1)-newOrigin-1;
        if syncshift ==-2,syncshift=0;end
        Data.load(syncDataPeriods,[],syncshift);
% $$$                 Data.load@MTAData(syncDataPeriods,[],syncshift);
% $$$                 fh = @(Data,syncDataPeriods,syncshift)load(Data,syncDataPeriods,[],syncshift);
% $$$                 feval(fh,Data,syncDataPeriods,syncshift)
        if ~isempty(find(syncZeroIndex,1)),
            Data.data(syncZeroIndex,:,:,:,:) = 0;
        end
    else
        if ~isempty(find(syncZeroIndex,1)),
            Data.data(syncZeroIndex,:,:,:,:) = 0;
        end
        nper = syncEpoch.copy;
        
        if (nper.size(1)-Data.size(1))==1, 
            nper.data = nper.data(1:Data.size(1));
            warning(['Synchronization offset 1 index larger than data']);
        end
        
        
        Data.data = Data.data(find(nper==1,1,'first'):find(nper==1,1,'last'),:,:,:,:);
    end
    
end

% Is this right? I think it is
if round(Data.origin*Data.sampleRate) ~= newOrigin;
    Data.origin = newOrigin/Data.sampleRate;
end

end
