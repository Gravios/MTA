function cax = apply_colorbar(axesHandle,location,colorscheme)
cax = colorbar(axesHandle,location);

cax.Units = 'centimeters';
colormap(axesHandle,colorscheme);

switch location
  case 'eastoutside'
    cax.Position(1) = cax.Position(1)+0.01;
    drawnow();
    pause(0.1);
    cax.Position(1) = sum(axesHandle.Position([1,3])+0.05);
  case 'southoutside'
    cax.Position = [sum(axesHandle.Position([1])),...
                    axesHandle.Position(2)-0.2,...
                    axesHandle.Position(3),...
                    0.15];
  otherwise 
    warning('MTA:utilities:graphics:apply_colorbar:UnknowLocation');
end
