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
Data = cat(1,flipud(Data(1:order,:,:,:,:,:)),...
                  Data,...
                  flipud(Data(end-order:end,:,:,:,:,:)));

dimensions = size(Data);
for n = 1:numApplications,
    Data(:,:,:,:,:,:) = ...
        reshape(...
         permute(...
          mean(circshift(GetSegs(Data(:,:,:,:,:,:),1:size(Data,1),order,nan),...
                         floor(order/2),2),...
               1,'omitnan'),...
                 [[2:numel(dimensions)+1],1]),...
                [size(Data,1),dimensions(2:end)]);
end

Data = Data(order+1:end-order-1,:,:,:,:,:);


% END MAIN -----------------------------------------------------------------------------------------