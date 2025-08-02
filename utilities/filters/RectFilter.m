function Data = RectFilter(Data,varargin)
% function Data = RectFilter(Data,varargin)
% 
% Inputs:
%  
%     Data:   (MTAData)      Data to be filtered (NxMxPx...) 
%     order:  (Integer)                   {default = 3} widow size in number of samples 
%     numApplications: (Integer)          number of passes with rectangular filter

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('order',                 3,                                                     ...
                 'numApplications',       3,                                                     ...
                 'dataType',              'linear'                                               ...
);
[order,numApplications,dataType] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------
nind = nniz(Data.data);

switch dataType
  case 'linear',
    Data.data = cat(1,flipud(Data.data(1:order,:,:,:,:,:)),...
                      Data.data,...
                      flipud(Data.data(end-order:end,:,:,:,:,:)));
  case 'circular'
    Data.data = cat(1,flipud(Data.data(end-order:end,:,:,:,:,:)),...
                      Data.data,...
                      flipud(Data.data(1:order,:,:,:,:,:)));
end

dimensions = size(Data);
for n = 1:numApplications,
    Data.data(nniz(Data),:,:,:,:,:) = ...
        reshape(...
         permute(...
          mean(circshift(GetSegs(Data(nniz(Data),:,:,:,:,:),1:sum(nniz(Data)),order,nan),...
                         -floor(order/2)-1,2),...
               1,'omitnan'),...
                 [[2:numel(dimensions)+1],1]),...
                [sum(nniz(Data)),dimensions(2:end)]);
end

Data.data = Data.data(order+1:end-order-1,:,:,:,:,:);
Data.data(~nind,:,:,:,:,:) = 0;

% UPDATE hash property of Data object
Data.update_hash(DataHash(struct('order',order,'numApplications',numApplications)));

% END MAIN -----------------------------------------------------------------------------------------