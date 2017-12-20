function Data = RectFilter(Data,varargin)
% function Data = RectFilter(Data,varargin)
% 
% Inputs:
%  
%     Data: (numeric:array,matrix)      Data to be filtered (NxMxPx...) 
%     order:  (Integer)                   {default = 3} widow size in number of samples 
%     numApplications: (Integer)          number of passes with rectangular filter

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('order',                 3,                                                     ...
                 'numApplications',       3                                                      ...
);
[order,numApplications] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------

dimensions = size(Data);
for n = 1:numApplications,
    Data.data(nniz(Data),:,:,:,:,:) = ...
        reshape(...
         permute(...
          mean(circshift(GetSegs(Data(nniz(Data),:,:,:,:,:),1:sum(nniz(Data)),order,nan),...
                         floor(order/2),2),...
               1,'omitnan'),...
                 [[2:numel(dimensions)+1],1]),...
                [sum(nniz(Data)),dimensions(2:end)]);
end

% UPDATE hash property of Data object
Data.update_hash(DataHash(struct('order',order,'numApplications',numApplications)));

% END MAIN -----------------------------------------------------------------------------------------