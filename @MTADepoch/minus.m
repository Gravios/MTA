function Data = minus(a,b)
% function Data = minus(a,b)
% 
% Add assertion that when two MTADepochs are inputed that
% they have the same origin

if isa(a,'MTADepoch')&&ismatrix(b)&&~isa(b,'MTADepoch'),
    Data = a.copy;
    if all(size(b)==1),
        Data.data = bsxfun(@minus,Data.data,b);
    else
        Data.data = SubstractRanges(Data.data,b);
    end
elseif isa(b,'MTADepoch')&&ismatrix(a)&&~isa(a,'MTADepoch'),
    Data = b.copy;
    if all(size(a)==1),
        Data.data = bsxfun(@minus,Data.data,a);
    else
        Data.data = SubstractRanges(a,Data.data);
    end
elseif isa(b,'MTADepoch')&&isa(a,'MTADepoch')
    b.resample(a.sampleRate);
    Data = a.copy;
    Data.data = SubstractRanges(Data.data,b.data);
    Data.label = [Data.label '-' b.label];
    Data.key = '';
end

end
