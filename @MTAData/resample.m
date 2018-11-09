function Data = resample(Data,DataObj,varargin)
% Data = resample(Data,DataObj)
% Resample Data to the target sample rate proveded in DataObj
% 
% NOTE - If target is an MTAData object rounding discrepancies will
%        be corrected by trimming the tails to fit the target object.
%    
% Input:
%   DataObject: numeric, target sample rate
%               MTAData, target object to match data's sampleRate
%
%   interpMethod: string, ['spline'] Method for data interpolation/ downsampling
%
% Assumes the two objects have their starting points synchronized
%
%    upsampling 
%         MTAData - uses spline
%         MTADepoch.type <- 'TimePeriods' - uses nearest neighbor
%
% WARNING - Doesn't modify the sync if sizes don't match
% WARNING - uses interp1 with spline as default
% 

%Data = Data.copy;
    
if isempty(DataObj), return; end
switch Data.type
  case 'TimeSeries'
    [interpMethod] = DefaultArgs(varargin,{'spline'},true);
    isMTA = isa(DataObj,'MTAData');
    if isMTA,
        if all(Data.size(1)==DataObj.size(1)) && Data.sampleRate==DataObj.sampleRate,
            return;
        end
        newSampleRate = DataObj.sampleRate;       
    else
        if Data.sampleRate==DataObj,
            return;
        end
        newSampleRate = DataObj;
    end

    Data.data(~nniz(Data),:,:,:,:) = 0;
    
    if isa(Data,'MTADepoch'),
        interpMethod = 'nearest';
    end
    
    if newSampleRate<Data.sampleRate&&~isa(Data,'MTADepoch'),
        Data.filter('ButFilter',3,newSampleRate/2.0000001,'low');
    end
    
    statusIsnumeric = isnumeric(Data.data);
    if ~statusIsnumeric, 
        dataClass = class(Data.data);
        Data.data = double(Data.data);
    end

% $$$     if newSampleRate<Data.sampleRate, shift = 0;%1./Data.sampleRate;
% $$$     else,                             shift = 0;
% $$$     end

    if isMTA
        ntvec = (1:DataObj.size(1))./newSampleRate;%+shift
        %ntvec = (1:DataObj.size(1))./newSampleRate;
        xtvec = (1:Data.size(1))./Data.sampleRate;
        Data.data = interp1(xtvec,Data.data,ntvec,interpMethod);
        Data.sampleRate = newSampleRate;
    else
        ntvec = (1:(Data.size(1)./Data.sampleRate.*newSampleRate))./newSampleRate;%+shift
        xtvec = (1:Data.size(1))./Data.sampleRate;
        Data.data = interp1(xtvec,Data.data,ntvec,interpMethod);
        Data.sampleRate = newSampleRate;
    end
    if Data.size(1)==1; Data.data = Data.data'; end

    if ~statusIsnumeric, 
        Data.data(~nniz(Data.data)) = 0;
        Data.data = feval(dataClass,Data.data);
    end
    
    if isa(Data,'MTADepoch'),
        Data.data(~nniz(Data.data)) = 0;
    end
    
    
  case 'TimePeriods'
    % Needs some more corrections for resampling

    if DataObj == 1
        newSampleRate = DataObj;
        rf = @(x)x;
        indshift = 0;
    elseif isa(DataObj,'MTAData')
        newSampleRate = DataObj.sampleRate;
        rf = @round;
    else
        newSampleRate = DataObj;
        rf = @round;
    end

    if newSampleRate == Data.sampleRate, return,end                
    
    if newSampleRate==1&&Data.sampleRate==1,
        indshift = 0;
    elseif  newSampleRate~=1&&Data.sampleRate==1,
        indshift = 1/newSampleRate;
    elseif  newSampleRate==1&&Data.sampleRate~=1,
        indshift = -1/Data.sampleRate;
    else
        indshift = 0;
    end                  


    Data.data = rf(Data.data/Data.sampleRate*newSampleRate+indshift);
    while sum(Data.data(:)==0)>1
        Data.data(1,:) = [];
    end
    
    if Data.size(1)>1,
        mind = find(Data.data(2:end,1)-Data.data(1:end-1,2)==0)+1;
    else
        mind = [];
    end
    
    if ~isempty(mind)
        Data.data(mind-1,2)=Data.data(mind,2);
        Data.data(mind,:)=[];
    end

    if isa(DataObj,'MTAData'),
        % Remove periods which exist past the end of the reference object
        while Data.data(end,1)>size(DataObj,1)
            Data.data(end,:) = [];
        end
        
        % Adjust final end point to match the end point of the reference object
        if Data.data(end)>size(DataObj,1)
            Data.data(end) = size(DataObj,1);
        end
    end

    if ~isempty(Data.data)&Data.data(1)==0,Data.data(1)=1;end

    
% $$$                   if isa(Data.sync,'MTAData'),
% $$$                       Data.origin = rf(Data.origin/Data.sampleRate*newSampleRate); 
% $$$                   else
% $$$                       Data.sync = rf(Data.sync/Data.sampleRate*newSampleRate);
% $$$                   end            
    Data.sampleRate = newSampleRate;
  otherwise
end
end

