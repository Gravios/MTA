function npoints = rotate_points_around_vectors(point, normalVector, angle)
%function mxyz = rotate_point_around_vector
%[angle,crossMarkers] = DefaultArgs(varargin,{45,{'hbx','hbr'}},1);
%
% NOTE: xyz must contain a marker named 'hcom' 
% NOTE: angle must be in degrees
% NOTE: refMarkers must only have two markers
%

nSamples = size(normalVector,1);

j =1:3;
head_kron = reshape(repmat(normalVector',3,1).*normalVector(:,j(ones(3,1),:)).',[3,3,nSamples]);
j = [ 0,-1, 1;...
      1, 0,-1;...
     -1, 1, 0];
k = [1,3,2;...
     3,1,1;...
     2,1,1];
head_cpm = reshape(normalVector(:,k)',3,3,nSamples) .* repmat(j,[1,1,nSamples]);

head_rotMat = cos(angle)*repmat(eye(3),[1,1,nSamples])...
              +sin(angle)*head_cpm...
              +(1-cos(angle))*head_kron;
j =1:3;
% Rotated marker;
npoints = permute(sum(head_rotMat.*permute(reshape(point(:,j(ones(3,1),:)),[nSamples,3,3]),[2,3,1]),2),[3,1,2]);




