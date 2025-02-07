
function newMarker = rotate_marker_around_vector(xyz,varargin)
%function mxyz = rotate_point_around_vector
%[angle,crossMarkers] = DefaultArgs(varargin,{pi/4,{'hbx','hbr'}},1);
%
% 
% NOTE: angle must be in radians

[angle,markers,oriMarker,refMarkers] = DefaultArgs(varargin,{pi,'head_right','hcom',{'head_back','head_right'}},1);

nind = nniz(xyz);

origin = xyz(nind,oriMarker,:);

points = sq(bsxfun(@minus,xyz(nind,markers,:), origin));

vec = bsxfun(@minus,xyz(nind,refMarkers,:), origin);
nvec = cross(sq(vec(:,1,:)),sq(vec(:,2,:)));
nvec = multiprod(nvec,1./sqrt(sum(nvec.^2,2)),2);

nmark = rotate_points_around_vectors(points, nvec, angle);

newMarker = xyz(:,markers,:);
newMarker(nind,:,:) = permute(nmark,[1,3,2]) + origin;
