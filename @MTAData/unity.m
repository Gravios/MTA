function [Data,meanA,stdA] = unity(Data,varargin)
%function [U,meanA,stdA] = unity(A,varargin)
%[ifnniz,meanA,stdA,drpOutPrctile,states] = DefaultArgs(varargin,{@nan,[],[],[],{}},true);
% A is n x m  matrix, and each column is a signal.
% the function gives a unity signal for each column
% i.e. [signal-mean(signal)] / std(signal)


% DEFARGS ------------------------------------------------------------------------------------------
defArgs = struct(...
    'ifnniz',            @nan,                ... (e.g. nan, zeros, ones)
    'meanA',             [],                  ...    
    'stdA',              [],                  ...
    'drpOutPrctile',     [],                  ...
    'states',            {{}},                ...
    'mean',              @nanmean,            ... (e.g. mean, nanmean, median)
    'std',               @nanstd,             ... (e.g. std, nanstd)       
    'profile',           '',                  ... unused  
    'overwrite',         false                ... unused  
);

[ifnniz,meanA,stdA,drpOutPrctile,states,...
 mean,std,profile,overwrite] = DefaultArgs(varargin,defArgs,'--struct');
%--------------------------------------------------------------------------------------------------




% MAIN ---------------------------------------------------------------------------------------------

nind = nniz(Data);

if ~isempty(states),
    states = MTADepoch.join(states);
    states.cast('TimeSeries',Data);
    states = logical(states.data);
else
    states = true([size(Data.data,1),1]);
end


A = Data.data;
if ~isempty(drpOutPrctile),
    newFeatureDomains = prctile(Data.data(nind,:,:,:,:,:),drpOutPrctile);
    inDomain = bsxfun(@lt,Data.data(:,:,:,:,:,:),newFeatureDomains(1,:)) & ...
               bsxfun(@gt,Data.data(:,:,:,:,:,:),newFeatureDomains(2,:));

    A(inDomain(:)) = nan;
end
A(A==0)=nan;
A(isinf(A))=nan;
A(~states,:) = nan;

if isempty(meanA),
    meanA = nanmean(A);
end
if isempty(stdA)
    stdA = nanstd(A);
end

Data.data = bsxfun(@rdivide,bsxfun(@minus,Data.data,meanA),stdA);

Data.data(~nind,:,:,:,:,:) = ifnniz([sum(~nind),size(Data,2)]);

% END MAIN -----------------------------------------------------------------------------------------