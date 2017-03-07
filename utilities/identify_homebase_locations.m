function [location,occupancy] = identify_homebase_locations(Trial)
%function homebaseLocation = identify_homebase_locations(Trial)
% Returns up to 3 locations which are to be homebases of rodent subject

xyz = preproc_xyz(Trial,'SPLINE_SPINE_HEAD_EQD');
state = Trial.stc{'pause'};
type = 'xy';
binSize = 20; %mm


[smoothedOccupancyMap,xyBins] = generate_occupancy_map(Trial,xyz,state,type,binSize);

normalizeOccupancyMap = smoothedOccupancyMap{1}./nansum(smoothedOccupancyMap{1}(:));

[location,occupancy] = LocalMinimaN(-normalizeOccupancyMap,...
                                    -max(normalizeOccupancyMap(:))/2,...
                                    200/diff(xyBins{1}(1:2)));
location(:,1) = xyBins{1}(location(:,1));
location(:,2) = xyBins{2}(location(:,2));

