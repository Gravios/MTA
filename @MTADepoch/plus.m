function Data = plus(a,b)
%function Data = plus(a,b)
%
% Inuput takes either a two element vector which indicates a time
%  shift to the start and stop point of each period within the
%  MTADepoch object.
%
% Note: if addition or subtraction results in a period where the
%       start precedes the stop then that period will be removed
%
% Note: if addition or subtraction results two periods overlaping
%       they will be merged into a single period
%
% Future: Addition of two MTADepochs will return their union
%       
% Future: suport TimeSeries versions

if isa(a,'MTADepoch')&&~isa(b,'MTADepoch')
    if strcmp(a.type,'TimePeriods'),
        if nargout==0,
            Data = a;
        else
            Data = a.copy;
        end

        if a.sampleRate>1,
            b = round(b*a.sampleRate);
        else
            b = b*a.sampleRate;
        end
        
        if prod(size(b) == Data.size)&&~(prod(size(b)==[1,2])),
            Data.data = Data.data+b;
        else
            Data.data = bsxfun(@plus,Data.data,b);
        end
        Data.clean();
    elseif strcmp(a.type,'TimeSeries'),
    end
    
elseif isa(b,'MTADepoch')&&~isa(a,'MTADepoch')
    if strcmp(a.type,'TimePeriods'),
        if nargout==0,
            Data = b;
        else
            Data = b.copy;
        end

        if b.sampleRate>1,
            a = round(a*b.sampleRate);        
        else 
            a = a*b.sampleRate;
        end
        
        if prod(size(a) == Data.size)&&~(prod(size(a)==[1,2])),
            Data.data = Data.data+a;
        else
            Data.data = bsxfun(@plus,Data.data,a);
        end
        
        Data.clean();
    elseif strcmp(a.type,'TimeSeries'),
    end
    
elseif isa(a,'MTADepoch')&&isa(b,'MTADepoch')
    if strcmp(a.type,'TimePeriods')&&strcmp(b.type,'TimePeriods'),
        if nargout==0,
            Data = a;
        else
            Data = a.copy;
        end
        Data.path = [];
        Data.filename = [];
        Data.label = [a.label '+' b.label];
        Data.key = 'c';
        if b.sampleRate ~= a.sampleRate, b.resample(a); end
        Data.data = JoinRanges(a.data,b.data);
        
        Data.clean();
    elseif strcmp(a.type,'TimeSeries'),
    end
    
end

% Join the new periods if they overlap
if sum((Data.data(2:end,1)-Data.data(1:end-1,2))<=0)>0;
    ndata = Data.data(1,:);
    for i = 2:Data.size(1),
        ndata = JoinRanges(ndata,Data.data(i,:));
    end
    Data.data = ndata;
end
if Data.sampleRate~=1,
    Data.data = round(Data.data);
end


end
