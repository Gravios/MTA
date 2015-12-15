gfunction [U,meanA,stdA] = unity(Data,varargin)
%function [U,meanA,stdA] = unity(A,varargin)
%[ifnniz,meanA,stdA,drpOutPrctile] = DefaultArgs(varargin,{@nan,[],[],[]},true);
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
        {}...
};

[ifnniz,meanA,stdA,drpOutPrctile,states] = DefaultArgs(varargin,defArgs,true);

nind = nniz(Data);

if ~isemtpy(states),
    states = join(states).cast('TimeSeries',Data);
    states = logical(states.data);
else
    states = ones(size(nind));
end

nind = nind&states;

if ~isempty(drpOutPrctile),
    Ao = prctile(Data(nind,:),drpOutPrctile);
else
    Ao = Data(nind,:);
end

if isempty(meanA)
    meanA = mean(Ao);
end
if isempty(stdA)
    stdA = std(Ao);
end

U = feval(ifnniz,size(A));

U(nind,:) = bsxfun(@rdivide,bsxfun(@minus,A(nind,:),meanA),stdA);


