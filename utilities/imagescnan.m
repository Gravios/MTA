function [handleImage, handleColorbar] = imagescnan(imageData, varargin)
% function [handleImage, handleColorbar] = imagescnan( imageData, varargin)
% imageData is either a Matrix to colormap plot or cell {Xaxis, Yaxis, Matrix}
% Alternative to imagesc with the ability to specify an hsv for NaN
% values
%
% Out:
%      handleImage:    image graphics handle
%      handleColorbar: colorbar graphics handle
%
% varargin:
%      colorLimits: default = min/max
%      dataType:                  'linear',  ...
%                 'colorbarIsRequired',        false,     ...
%                 'nanRGB',                    [1,0,0.8], ...
%                 'gamma',                     1,         ...
%                 'value',                     1,         ...
%                 'colorMap',                  @parula    ...
%
% if 'sym' then will make it +/- max(min,max)
% nanRGB default = [1 0 0.8] (light gray)
% dataIsCircular: a flag that defines whether the color map should be circular or
% not 
%
% Value - scales the hue according to the value
%
% Evgeny Resnik: Additional check if colorLimits is empty when
%                imageData has only NaN values
% Evgeny Resnik: Added h2 as an optional output parameters to be able to
%                change properties of colorbars (06.2012).
% DEFARGS ------------------------------------------------------------------------------------------
defargs = struct('colorLimits',               [],        ...
                 'dataType',                  'linear',  ...
                 'colorbarIsRequired',        false,     ...
                 'nanRGB',                    [0.5,0.5,0.5], ...
                 'gamma',                     1,         ...
                 'value',                     1,         ...
                 'colorMap',                  @parula    ...
);
[colorLimits,dataType,colorbarIsRequired,nanRGB,...
 gamma,value,colorMap] = DefaultArgs(varargin,defargs,'--struct');
%--------------------------------------------------------------------------------------------------

% MAIN ---------------------------------------------------------------------------------------------

if iscell(imageData)
% UNPACK cell into variables
    Xax = imageData{1};
    Yax = imageData{2};
    imageData = imageData{3};
else
% ASSIGN image axes' domains
    Xax = [1:size(imageData,1)];
    Yax = [1:size(imageData,2)];
end

if ~isreal(imageData)
    Value = abs(imageData);
    Value = (Value-min(Value(:)))/(max(Value(:))-min(Value(:)));
    imageData = angle(imageData);
    dataType = 'complex';
end

if isempty(colorLimits)
% SET color scale domain
   colorLimits =  [min(min(imageData(~isnan(imageData)))),...
                   max(max(imageData(~isnan(imageData))))];
elseif strcmp(colorLimits,'sym'),
    colorLimits = repmat(max(abs([min(min(imageData(~isnan(imageData)))),...
                                  max(max(imageData(~isnan(imageData))))])),...
                         [1,2])*0.7;
end

if isempty(colorLimits),  colorLimits = [0 1]; end

switch dataType
  case 'linear'
    cm = colorMap(1000);
    bins = discretize(imageData,linspace([colorLimits,1000]));
    finalImage = reshape(repmat(bins,[1,1,3]),[],3);
    finalImage(nniz(finalImage),:) = cm(finalImage(nniz(finalImage),1),:);
    finalImage(isnan(imageData(:)),:) = repmat(nanRGB,[sum(isnan(imageData(:))),1]);
    finalImage = reshape(finalImage,[numel(Yax),numel(Xax),3]);


  case 'circular'
    cm = colorMap(1000);
    bins = discretize(imageData,linspace([colorLimits,1000]));
    finalImage = reshape(repmat(bins,[1,1,3]),[],3);
    finalImage(nniz(finalImage),:) = cm(finalImage(nniz(finalImage),1),:);
    finalImage(isnan(imageData(:)),:) = repmat(nanRGB,[sum(isnan(imageData(:))),1]);
    finalImage = reshape(finalImage,[numel(Yax),numel(Xax),3]);
  
  case 'complex', ... may need some work
    Hsv(:,:,1) = clip((imageData-colorLimits(1))./(colorLimits(2)-colorLimits(1)),0,1).^Gamma;
    % for complex data
    % Amp = abs(C);
    % AmpMin = min(Amp(:));
    % AmpMax = max(Amp(:));
    % Phase = angle(C);
    % 
    % Hsv = zeros([size(C) 1]);
    % 
    % Hsv(:,:,1) = mod(Phase/(2*pi), 1);
    % Hsv(:,:,2) = 1;
    % 
    % if abs(AmpMax - AmpMin) < 1e-5 * max(abs(AmpMax), abs(AmpMin));
    % 	Hsv(:,:,3) = 1;
    % else
    % 	Hsv(:,:,3) = ((Amp-AmpMin)/(AmpMax-AmpMin)).^gamma;
    % end
    
end



try 
% DISPLAY image in current axis
    image(Xax,Yax,finalImage);
catch
    fprintf('\nimage(hsv2rgb(Hsv)) failed. Trying image(abs(hsv2rgb(Hsv)))\n')
    image(Xax,Yax,abs(finalImage));
end

set(gca, 'ydir', 'reverse')

handleImage    = gca;
handleColorbar = SideBar;
image(0,linspace([colorLimits,1000]), permute(cm, [1,3,2]));
set(handleColorbar, 'ydir', 'normal');
set(handleColorbar, 'xtick', []);
set(handleColorbar, 'yaxislocation', 'right');
axes(handleImage);


% END MAIN -----------------------------------------------------------------------------------------