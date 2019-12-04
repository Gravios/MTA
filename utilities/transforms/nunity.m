%function [U,meanA,stdA] = nunity(A,varargin)
%[ifnniz,meanA,stdA,drpOutPrctile] = DefaultArgs(varargin,{@nan,[],[],[]},true);
% A is n x m  matrix, and each column is a signal.
% the function gives a unity signal for each column
% i.e. [signal-mean(signal)] / std(signal)

function [U,meanA,stdA] = nunity(A,varargin)
[ifnniz,meanA,stdA,drpOutPrctile,refDim] = DefaultArgs(varargin,{@nan,[],[],[],[]},true);

nind = nniz(A);

if ~isempty(drpOutPrctile),
    boundaries = prctile(A(nind,:,:,:,:,:),drpOutPrctile);
    Ao = reshape(A(repmat(nind,[1,size(A,2),size(A,3),size(A,4),size(A,5),size(A,6)]) ...
           & A > repmat(boundaries(1,:,:,:,:,:),[size(A,1),1,1,1,1,1]) ...
           & A < repmat(boundaries(2,:,:,:,:,:),[size(A,1),1,1,1,1,1])),...
                 [],size(A,2),size(A,3),size(A,4),size(A,5),size(A,6));
else
    Ao = A(nind,:,:,:,:,:);
end


if isempty(meanA)
    if ~isempty(refDim),
        meanA = repmat(mean(Ao(:,refDim,:,:,:)),[1,size(Ao,2),size(Ao,3),size(Ao,4),size(Ao,5)]);
    else
        meanA = mean(Ao);
    end    
end
if isempty(stdA)
    if ~isempty(refDim),
        stdA = repmat(std(Ao(:,refDim,:,:,:)),[1,size(Ao,2),size(Ao,3),size(Ao,4),size(Ao,5)]);
    else
        stdA = std(Ao);
    end
end


U = feval(ifnniz,size(A));

U(nind,:,:,:,:,:) = bsxfun(@rdivide,bsxfun(@minus,A(nind,:,:,:,:,:),meanA),stdA);

