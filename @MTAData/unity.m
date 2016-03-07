function [Data,meanA,stdA] = unity(Data,varargin)
%function [U,meanA,stdA] = unity(A,varargin)
%[ifnniz,meanA,stdA,drpOutPrctile,states] = DefaultArgs(varargin,{@nan,[],[],[],{}},true);
% A is n x m  matrix, and each column is a signal.
% the function gives a unity signal for each column
% i.e. [signal-mean(signal)] / std(signal)


defArgs = {...
 ...    ifnniz, (e.g. @nan, @zeros, @ones)
        @nan,          ...
 ...
 ...    meanA    
        [],            ...
 ...
 ...    stdA
        [],            ...
 ...
 ...    drpOutPrctile
        [],            ...
 ...
 ...    states
        {},            ...
 ...
 ...    mean (e.g. mean, nanmean, median)
        @nanmean,      ...
 ...
 ...    std  (e.g. std, nanstd)       
        @nanstd        ...
};

[ifnniz,meanA,stdA,drpOutPrctile,states,mean,std] = DefaultArgs(varargin,defArgs,true);

nind = nniz(Data);

if ~isempty(states),
    states = MTADepoch.join(states);
    states.cast('TimeSeries',Data);
    states = logical(states.data);
else
    states = true(size(nind));
end

nind = nind&states;

if ~isempty(drpOutPrctile),
    A = prctile(Data.data(nind,:),drpOutPrctile);
else
    A = Data.data(nind,:);
end

if isempty(meanA)
    meanA = mean(A);
end
if isempty(stdA)
    stdA = std(A);
end

Data.data = feval(ifnniz,size(A));

Data.data(nind,:) = bsxfun(@rdivide,bsxfun(@minus,A,meanA),stdA);


