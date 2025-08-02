function rcsd = csd(Data,varargin)

% DEFARGS ------------------------------------------------------------------------------------------

defargs = struct('channelInterval' , 1,                                                          ...
                    'channelPitch' , 50 );

[channelInterval,channelPitch] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------

numCsdChannels = size(Data.data,2)-2.*channelInterval-1;
rcsd = copy(Data);
rcsd.data = Data.data(:,1:numCsdChannels);
for c = 1:numCsdChannels;
    rcsd.data(:,c) = (Data.data(:,c)+Data.data(:,c+3)-2.*Data.data(:,c+1))./(2.*channelInterval.*channelPitch);
end

