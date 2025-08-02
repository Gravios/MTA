% toggle_axes_deletion

global DELETE_CURRENT_AXES

DELETE_CURRENT_AXES = ~DELETE_CURRENT_AXES;

if DELETE_CURRENT_AXES
    warning(['MTA:utilities:graphics:toggle_axes_deletion:DELETE_CURRENT_AXES is True']);
else
    warning(['MTA:utilities:graphics:toggle_axes_deletion:DELETE_CURRENT_AXES is False']);
end
