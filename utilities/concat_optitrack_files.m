function [xyzData,markerData,xyzSampleRate] = concat_optitrack_files(Session)
%[xyzData,markers] = concat(Session)
% Determine number of position tracking 
% trials and concantenates parts
dirArray = dir(fullfile(Session.spath, Session.maze.name));
trialArray = {};
trialPartSize = {};
lastTrial = 0;
expr = 'trial..._';
for i = 1:length(dirArray);
    validTrial = regexpi(dirArray(i).name,expr);
    if validTrial,
        validTrial = str2num(dirArray(i).name([validTrial+5:validTrial+7]));
        if validTrial>lastTrial,
            lastTrial = validTrial;
            trialArray{lastTrial,1} = dirArray(i).name;
            trialPartSize{lastTrial} = 1;
        else
            trialPartSize{lastTrial} = trialPartSize{lastTrial} + 1;
            trialArray{lastTrial,trialPartSize{lastTrial}} = dirArray(i).name;
        end
    end
end
% Sort the names within the cells
% This is assumed to be in order for the moment
%Load Files and concantenate
xyzData = cell(lastTrial,max([trialPartSize{:}]));
xyzSampleRate = [];
number_of_markers = 0;
for j = 1:lastTrial,
    for k = 1:trialPartSize{j},
        load(fullfile(Session.spath, Session.maze.name, trialArray{j,k}));
        % Check if all parts have the same number of
        % markers as the first part
        if length(markers)~=number_of_markers&&(j~=1&&k~=1),
            warning(['Number of markers changes at part: ' trialArray{j,k}])
        else
            number_of_markers = length(markers);
        end
        if exist('sampleRate','var'),
            xyzSampleRate(end+1) = sampleRate; 
        end
        xyzData{j,k} = xyzpos;
        markerData{j,k} = markers;
    end
end
xyzSampleRate = unique(xyzSampleRate);
%assert(numel(xyzSampleRate)<=1,'MTA:utilities:concatViconFiles: Sample rate changes between files/fileparts');
end
