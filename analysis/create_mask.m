function [mask] = create_mask(bins,shape)

    width = numel(bins{1});    
    height = numel(bins{2});

    
    radius = round(numel(pfsBins{1})/2)-find(pfsBins{1}<-440,1,'last');
    centerW = width/2;
    centerH = height/2;
    [W,H] = meshgrid(1:width,1:height);           
    circMask = repmat(double(sqrt((W-centerW-.5).^2 + (H-centerH-.5).^2) < radius),[1,1,pfsBinsDims(3:4)]);
