function imd = interMarkerDistance(xyz)
%function imd = interMarkerDistance(xyz)
%finds the pairwise distances between markers
%
% Input:
%     xyz: MTADxyz, Timeseries containing the xyz dimensions of 
%                   a set of marker in space

imd = sqrt(sum( markerDiffMatrix(xyz).^2,4));
