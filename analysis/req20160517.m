function req20160517(Trial,varargin)
Trial = 'jg05-20120317.cof.all';
Trial = MTATrial.validate(Trial);
xyz = Trial.load('xyz');



markers = {'head_back','head_left','head_front','head_right'};
A = xyz(:,markers,:);
A = bsxfun(@minus,A,mean(A));
A = permute(A,[2,3,1]);
A = sq(mat2cell(A,size(A,1),size(A,2),ones([size(A,3),1])));
[U,S,V]=cellfun(@svd,A,'uniformoutput',false);

[t,p,r] = cart2sph(coeff(1),coeff(2),coeff(3))



V1=[coeff(1:3)]';

V2=[0 -1 0];
rotationAxis=cross(V1,V2);
rotationAxis=rotationAxis/norm(rotationAxis);
 
%Calculate the theta
theta=acos(dot(V1,V2)/(norm(V1)*norm(V2)));

%R= rotationmat3D(theta,rotationAxis);

figure,plotSkeleton(Trial,xyz,t);