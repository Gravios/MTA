function diffMat = markerDiffMatrix(xyz)
%diffMat = markerDiffMatrix(Session)
%Create a time series where the position of every marker
%has been substracted from one another
%
%  Input
%
%  Output: 
%
%      diffMat: numericArray, (index,marker1,marker2,dim) Each markers
%                              position is substracted from every others
%                              to create all marker based frames of reference.

nframe = size(xyz,1); %nframe: number of frames (time)
nmar   = size(xyz,2);  %nmar: number markers 
ndim   = size(xyz,3); %ndim: number of spatial dimensions (xyz)
j =1:nmar;
diffMat = permute(cat(4,permute(reshape(repmat(xyz(:,:,1)',nmar,1)-xyz(:,j(ones(nmar,1),:),1).',[nmar,nmar,nframe]),[3,1,2]),permute(reshape(repmat(xyz(:,:,2)',nmar,1)-xyz(:,j(ones(nmar,1),:),2).',[nmar,nmar,nframe]),[3,1,2]),permute(reshape(repmat(xyz(:,:,3)',nmar,1)-xyz(:,j(ones(nmar,1),:),3).',[nmar,nmar,nframe]),[3,1,2])),[1,3,2,4]);

