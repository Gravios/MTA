% function U = unity(A) 
% A is n x m  matrix, and each column is a signal.
% the function gives a unity signal for each column
% i.e. [signal-mean(signal)] / std(signal)

function [U,meanA,stdA] = nunity(A,varargin)
[ifnniz,meanA,stdA] = DefaultArgs(varargin,{@nan,[],[]},true);

nind = nniz(A);
if isempty(meanA)
    meanA = mean(A(nind,:));
end
if isempty(stdA)
    stdA = std(A(nind,:));
end


U = feval(ifnniz,size(A));

U(nind,:) = bsxfun(@rdivide,bsxfun(@minus,A(nind,:),meanA),stdA);

