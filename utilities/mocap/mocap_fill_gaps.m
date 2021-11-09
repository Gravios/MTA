function obj = mocap_fill_gaps(obj,gapPeriods,minGapLength)

gapLen = diff(gapPer,1,2);

gapConv = zeros([size(obj,1),1]);
gapConvPer = gapPeriods(ismember(gapLen,1:minGapLen),:);
for gid = 1:size(gapConvPer,1)
gapConv(gapConvPer(gid,1):gapConvPer(gid,2)) = 1;
end
gapConv = conv(gapConv,ones([minGapLen,1]),'same');

gapPer = ThreshCross(gapConv,0.1,0);
gapLen = diff(gapPer,1,2);
gapPer = bsxfun(@plus,gapPer,[1,0]);
gapDstBA = [gapPer(:,1);size(obj,1)]-[1;gapPer(:,2)];

for gid = find(gapDstBA(1:end-1)>2.*gapLen & gapDstBA(2:end)>2.*gapLen)'
    
    halfSegLen = round(sum([1,2].*gapLen(gid)));
    segInd = gapPer(gid,1)-halfSegLen:gapPer(gid,2)+halfSegLen;
    seg = obj(segInd,1,1);
    nanVec  = ~nniz(seg);
    posInds = find(~nanVec);
    nanInds = find(nanVec); 
    
    for rid = 1:size(obj,2)
        for did = 1:size(obj,3)
            obj.data(segInd(nanInds),rid,did) = interp1(posInds,obj(segInd(posInds),rid,did),nanInds,'pchip');
        end
    end
end
