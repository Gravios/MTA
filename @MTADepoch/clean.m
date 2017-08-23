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
Data.data = sort(Data.data);

% Is this if statement really necessary  
if Data.sampleRate==1,
% DROP periods which exist before or conains the origin
    Data.data(Data.data(:,1)<0|Data.data(:,2)<0,:) = [];
% DROP periods which exceed the end of the sync
    Data.data(Data.data(:,1)>Data.sync.data(end)) = [];
    if Data.isempty(), return, end    
% TRUNCATE periods which terminate after end of the sync
    Data.data( Data.data(:,1)<Data.sync.data(end)...
        &Data.data(:,2)>Data.sync.data(end),2)...
        = Data.sync.data(end);        
% MERGE and remove periods which are overlapping or joint    
    Data.data = [Data.data(:,1),circshift(Data.data(:,2),1)];
    mpss = Data.data(1,:);
    Data.data(-diff(Data.data,1,2)<=0,:)=[];
    Data.data = cat(1,mpss,Data.data);
    Data.data = [Data.data(:,1),circshift(Data.data(:,2),-1)];
    
else
% DROP periods which exist before or conains the origin
    Data.data(Data.data(:,1)<=0|Data.data(:,2)<=0,:) = [];
% DROP periods which exceed the end of the sync
    Data.data(Data.data(:,1)>round(Data.sync.data(end)*Data.sampleRate)) = [];
% TRUNCATE periods which terminate after end of the sync
    if Data.isempty(), return, end
    Data.data( Data.data(:,1)<round(Data.sync.data(end)*Data.sampleRate)...
        &Data.data(:,2)>round(Data.sync.data(end)*Data.sampleRate),2)...
        = round(Data.sync.data(end)*Data.sampleRate);
% MERGE and REMOVE periods which are overlapping or joint
    Data.data = [Data.data(:,1),circshift(Data.data(:,2),1)];
    mpss = Data.data(1,:);
    Data.data(-diff(Data.data,1,2)<=0,:)=[];
    Data.data = cat(1,mpss,Data.data);
    Data.data = [Data.data(:,1),circshift(Data.data(:,2),-1)];

end

end
