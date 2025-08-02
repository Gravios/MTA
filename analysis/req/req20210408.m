% req20210408
%      Tags: quaternion
%      Status: Active
%      Type: Utility
%      Author: Justin Graboski
%      Final_Forms: NA
%      Project: General
%      Description: building out quaternion functions

% Quaternion time
%    1   i   j   k
% 1  1   i   j   k
% i  i  -1   k  -j
% j  j  -k  -1   i
% k  k   j  -i  -1
% 
% ijk = -1
%
%sum( [1,6,11,16] .* [1,-1,-1,-1] )
%sum( [2,5,12,15] .* [1, 1,-1, 1] )
%sum( [3,8, 9,14] .* [1, 1, 1,-1] )
%sum( [4,7,10,13] .* [1,-1 ,1, 1] )

hpi = [1, 6,11,16, 2, 5,12,15, 3, 8, 9,14, 4, 7,10,13];
hps = [1,-1,-1,-1, 1, 1,-1, 1, 1, 1, 1,-1, 1,-1 ,1, 1];

ind = [Trial.stc{'gper'}];
quat = sq(Arena(ind,'Arena',5:8));
quatC = bsxfun(@times,quat,[1,-1,-1,-1]);
quatR = cat(2,zeros([size(quat,1),1]),sq(ratAC(ind,1,1:3)));

tic
quatM = multiprod(permute(quatR,[1,3,2]),quatC,[2,3],[2,3]);
quatM = multiprod(quatR,quatC,[2,0],[0,2]);
quatM = reshape(quatM,[size(quatM,1),16]);
quatM = permute(sum(reshape(bsxfun(@times,quatM(:,hpi),hps),[size(quatM,1),4,4]),2),[1,3,2]);
quatM = multiprod(quat,quatM,[2,0],[0,2]);
quatM = reshape(quatM,[size(quatM,1),16]);
quatM = permute(sum(reshape(bsxfun(@times,quatM(:,hpi),hps),[size(quatM,1),4,4]),2),[1,3,2]);
toc

quatM = bsxfun(@times,quatM,[0,-1, 1,-1]);

figure();
hold('on');
plot3(quatR(:,2),quatR(:,3),quatR(:,4),'.b')
plot3(quatM(:,2),quatM(:,3),quatM(:,4),'.r')