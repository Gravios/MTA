function R = SubstractRanges(R1,R2)
%function R = SubstractRanges(R1,R2)
%
% substracts from ranges R2 ranges from R1 : R = R1-R2

NotR2 = [[-inf, R2(1,1)];[R2(1:end-1,2), R2(2:end,1)]; [R2(end,2), inf]];
R = IntersectRanges(R1, NotR2);

