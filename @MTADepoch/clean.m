function Data = clean(Data)
% function Data = clean(Data)
% Internal utility to clean up periods after set operations
%
% Drop periods with zero or negative duration. eg([5,3] or [[10,10])
% Drop periods which exist before or conains the origin
% Drop periods which exceed the end of the sync
% Truncate periods which terminate after end of the sync
%
% Drop periods with zero or negative duration. eg([5,3] or [[10,10])
perDur = diff(Data.data,1,2);
Data.data(perDur<=0,:) = [];
% Drop periods which exist before or conains the origin
Data.data(Data.data(:,1)<=0|Data.data(:,2)<=0,:) = [];
% Drop periods which exceed the end of the sync
Data.data(Data.data(:,1)>round(Data.sync.data(end)*Data.sampleRate)) = [];
% Truncate periods which terminate after end of the sync
Data.data( Data.data(:,1)<round(Data.sync.data(end)*Data.sampleRate)...
           &Data.data(:,2)>round(Data.sync.data(end)*Data.sampleRate),2)...
    = round(Data.sync.data(end)*Data.sampleRate);
end
