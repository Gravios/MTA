function rotmat = transform_vector_to_rotation_matrix(xyz, markerPair, rotation)
%function rotmat = transform_vector_to_rotation_matrix(xyz, markerPair, rotation)
%
% Create transformation matrix for a vector pair from a MTAData
% object.

rotmat = xyz(:,markerPair{2},[1,2])-xyz(:,markerPair{1},[1,2]);
rotmat = sq(bsxfun(@rdivide,rotmat,sqrt(sum(rotmat.^2,3))));
rotmat = cat(3,rotmat,sq(rotmat)*[0,-1;1,0]);
rotmat = multiprod(rotmat,...
                   [cos(rotation),-sin(rotation); ...
                    sin(rotation), cos(rotation)],...
                   [2,3],...
                   [1,2]);
