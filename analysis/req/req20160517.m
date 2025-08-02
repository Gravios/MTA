function req20160517(Trial,varargin)
Trial = 'jg05-20120317.cof.all';
%Trial = 'er02-20110908.cof.all';
Trial = MTATrial.validate(Trial);
xyz = Trial.load('xyz');



markers = {'head_back','head_left','head_front','head_right'};
A = xyz(:,markers,:);
A = bsxfun(@minus,A,mean(A));
A = permute(A,[2,3,1]);
A = sq(mat2cell(A,size(A,1),size(A,2),ones([size(A,3),1])));
[U,S,V]=cellfun(@svd,A,'uniformoutput',false);

t = zeros([numel(V),1]);
p = zeros([numel(V),1]);
r = zeros([numel(V),1]);

mPlane = zeros([numel(V),3]);

for i = 1:numel(V),
    mPlane(i,:) = V{i}(:,3);
end

mPlane(mPlane(:,3)<0,:) = sign(mPlane(mPlane(:,3)<0,:)).*mPlane(mPlane(:,3)<0,:);

for i = 1:numel(V),
    [t(i),p(i),r(i)] = cart2sph(mPlane(i,1),mPlane(i,2),mPlane(i,3));
end


figure,plot(asin(sqrt(sum(cross(mPlane,circshift(mPlane,-1)).^2,2))))

th = t;
th(p<0)=th(p<0)+pi/2;
th(th>pi) = th(th>pi)-pi;

V1=[coeff(1:3)]';

V2=[0 -1 0];
rotationAxis=cross(V1,V2);
rotationAxis=rotationAxis/norm(rotationAxis);
 
%Calculate the theta
theta=acos(dot(V1,V2)/(norm(V1)*norm(V2)));

%R= rotationmat3D(theta,rotationAxis);

figure,plotSkeleton(Trial,xyz,t);