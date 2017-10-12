function signal = RectFilter(signal,varargin)
% function signal = RectFilter(signal,varargin)
% 
% Inputs:
%  
%     signal: (numeric:array,matrix)      Signal to be filtered (NxMxPx...) 
%     order:  (Integer)                   {default = 3} widow size in number of samples 
%     numApplications: (Integer)          number of passes with rectangular filter

% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('order',                 3,                                                     ...
                 'numApplications',       3                                                      ...
);
[order,numApplications] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------



% MAIN ---------------------------------------------------------------------------------------------

signalDimensions = size(signal);
for n = 1:numApplications,
    signal(nniz(signal),:,:,:,:,:) = ...
        reshape(...
         permute(...
          mean(circshift(GetSegs(signal(nniz(signal),:,:,:,:,:),1:sum(nniz(signal)),order,nan),...
                         floor(order/2),2),...
               1,'omitnan'),...
                 [[2:ndims(signal)+1],1]),...
                [sum(nniz(signal)),signalDimensions(2:end)]);
end

% END MAIN -----------------------------------------------------------------------------------------