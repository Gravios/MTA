function Data = fillgaps(Data,varargin)
% function Data = fillgaps(Data,gap_size)
% fills gaps between epochs which are smaller than the specified
% size
%
% varargin:
%   gap_size: numeric, gap size in seconds
%
% NOTE: Only functions for TimePeriods
%
if ~isempty(varargin),
    gap_size = varargin{1};
else
    gap_size = round(.1*Data.sampleRate);
end


switch Data.type
  case 'TimePeriods'
    if size(Data.data,1)>1,
        interPerDur = Data.data(2:end,1)-Data.data(1:end-1,2);
        c = 1;
        while ~isempty(interPerDur)
            if interPerDur(1)<gap_size,
                Data.data(c,:)   = [Data.data(c,1),Data.data(c+1,2)];
                Data.data(c+1,:) = [];
            else
                c = c+1;
            end
            interPerDur(1) = [];
        end
    end
  case 'TimeSeries'
    perStart = find(diff(double(Data.data))== 1)+1;
    perStop =  find(diff(double(Data.data))== -1);
    
    if ~isempty(perStart)&&~isempty(perStop),
        if perStop(1)<perStart(1),perStart = [1;perStart]; end
        if perStop(end)<perStart(end),perStop = [perStop;size(Data,1)]; end
        
        interPerDur = perStart(2:end)-perStop(1:end-1);
        c = 1;
        while ~isempty(interPerDur)
            if interPerDur(1)<gap_size,
                Data.data(perStop(c):perStart(c+1)) = 1;
                Data.data(c+1,:) = [];
            else
                c = c+1;
            end
            interPerDur(1) = [];
        end
    end
  otherwise
    error('MTA:MTADepoch:WhatDidYouDO!!!')

end
end


