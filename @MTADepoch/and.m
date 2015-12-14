function Data = and(a,b)
%function Data = and(a,b)
%
%If one of the imputs is a normal array then
%the array elements are assumed to be in the 
%sampling rate of the MTADepoch object.
%
% TODO - add MTADepoch/MTADepoch comparisions
Data = [];
if isa(a,'MTADepoch')&&ismatrix(b)
    Data = a.copy;
    switch a.type
      case 'TimePeriods'
        Data.data = IntersectRanges(Data.data,b);
        perDur = diff(Data.data,1,2);
        Data.data(perDur<=0,:) = [];
      case 'TimeSeries'
        Data.data = Data.data&b;
    end
    
elseif isa(b,'MTADepoch')&&ismatrix(a)
    Data = b.copy;
    switch b.type
      case 'TimePeriods'
        Data.data = IntersectRanges(Data.data,a);
        perDur = diff(Data.data,1,2);
        Data.data(perDur<=0,:) = [];
      case 'TimeSeries'
        Data.data = Data.data&a;
    end

elseif isa(b,'MTADepoch')&&isa(a,'MTADepoch')
    %Data = 
end
end
