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
                 'numApplications',       3,                                                     ...
                 'dataType',              'linear'                                               ...
);
[order,numApplications,dataType] = DefaultArgs(varargin,defargs,'--struct');
%---------------------------------------------------------------------------------------------------


% MAIN ---------------------------------------------------------------------------------------------
dimensions = size(Data);

% $$$ figure,plot(sq(Data(:,20,101,1,1)))
% $$$ tData = Data;

for n = 1:numApplications,
    switch dataType
      case 'linear',
        Data = cat(1,flipud(Data(1:order,:,:,:,:,:)),...
                   Data,...
                   flipud(Data(end-order+1:end,:,:,:,:,:)));
      case 'circular'
        Data = cat(1,Data(end-order+1:end,:,:,:,:,:),...
                   Data,...
                   Data(1:order,:,:,:,:,:));
      otherwise
        error('MTA:utilities:@double:RectFilter:UnknownDataType');
    end

    %hold('on');plot(-3:19,sq(Data(:,20,101,1,1)))

    Data = circshift(permute(mean(GetSegs(Data,1:size(Data,1),order,nan),1,'omitnan'),[2:numel(dimensions)+1,1]),-floor(order/2)-1);
    Data = Data(1:dimensions(1),:,:);
                    

% $$$     Data(:,:,:,:,:,:) = ...
% $$$         reshape(...
% $$$          permute(...
% $$$           mean(circshift(GetSegs(Data(:,:,:,:,:,:),1:size(Data,1),order,nan),...
% $$$                          floor(order/2)+1,2),...
% $$$                1,'omitnan'),...
% $$$                  [[2:numel(dimensions)+1],1]),...
% $$$                 [size(Data,1),dimensions(2:end)]);
end
Data = reshape(Data,dimensions);
%hold('on');plot(sq(Data(:,20,101,1,1)))


% END MAIN -----------------------------------------------------------------------------------------