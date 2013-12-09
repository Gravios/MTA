%function R = SubstractRanges(R1,R2)
% substracts from ranges R1 ranges from R2 : R = R1-R2

function R = SubstractRanges(R1,R2)


NotR2 = [R2(1:end-1,2) R2(2:end,1); [R2(end,2) inf]];
if R2(1,1)>1
    NotR2 =[[1 R2(1,1)]; NotR2];
end

R = IntersectRanges(R1, NotR2);

