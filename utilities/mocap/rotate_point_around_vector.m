function mxyz = rotate_point_around_vector(xyz,varargin)
%function mxyz = genRotatedMarker(xyz,varargin)
%[angle,crossMarkers] = DefaultArgs(varargin,{45,{'hbx','hbr'}},1);
%
% NOTE: xyz must contain a marker named 'hcom' 
% NOTE: angle must be in degrees
% NOTE: refMarkers must only have two markers
%

[markers,angle,refMarkers] = DefaultArgs(varargin,{'hbx',45,{'hbx','hrx'}},1);

nind = nniz(xyz);
mxyz = xyz(:,'hcom',:)
xyz_t = sq(xyz(nind,markers,:)-mxyz);

xyz_hb = bsxfun(@minus,xyz(:,refMarkers,:),xyz(:,'hcom',:));
head_norm = cross(sq(xyz_hb(:,1,:)),sq(xyz_hb(:,2,:)));
head_norm = multiprod(head_norm,1./sqrt(sum(head_norm.^2,2)),2);
j =1:3;

head_kron = reshape(repmat(head_norm',3,1).*head_norm(:,j(ones(3,1),:)).',[3,3,size(head_norm,1)]);
j = [ 0,-1, 1;...
      1, 0,-1;...
     -1, 1, 0];
k = [1,3,2;...
     3,1,1;...
     2,1,1];
head_cpm = reshape(head_norm(:,k)',3,3,size(head_norm,1))...
           .*repmat(j,[1,1,size(head_norm,1)]);
rot_ang = deg2rad(angle);
head_rotMat = cos(rot_ang)*repmat(eye(3),[1,1,size(head_norm,1)])...
              +sin(rot_ang)*head_cpm...
              +(1-cos(rot_ang))*head_kron;
j =1:3;
% Rotated marker;
nmark = permute(sum(head_rotMat.*permute(reshape(xyz_t(:,j(ones(3,1),:)),[size(head_norm,1),3,3]),[2,3,1]),2),[3,1,2]);


mxyz(nind,1,:) = permute(nmark,[1,3,2]) + mxyz;

