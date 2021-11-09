function Data = fill_gaps(Data,minGapLength,maxGapLength)
% function Data = fill_gaps(Data,minGapLength,maxGapLength)
% 
% Someone yell at me for not documenting this
% 

nanPer = ThreshCross(isnan(Data.data(:,1,1)),0.5,0);
gapLen = diff(nanPer,1,2);
nanPer = bsxfun(@plus,nanPer,[1,0]);
gapDstBA = [nanPer(:,1);size(Data.data,1)]-[1;nanPer(:,2)];
% find smallest gaps with no gaps to the before and after by distance 2*length(gap)

for gid = find(gapDstBA(1:end-1)>2.*gapLen & gapDstBA(2:end)>2.*gapLen)'
    
    halfSegLen = round(sum([1,2].*gapLen(gid)));
    segInd = nanPer(gid,1)-halfSegLen:nanPer(gid,2)+halfSegLen;
    seg = Data.data(segInd,1,1);
    nanVec  = ~nniz(seg);
    posInds = find(~nanVec);
    nanInds = find(nanVec); 
    
    for rid = 1:size(Data.data,2)
        for did = 1:size(Data.data,3)
            Data.data(segInd(nanInds),rid,did) = interp1(posInds,Data.data(segInd(posInds),rid,did),nanInds,'linear');
        end
    end
end

gapPer = ThreshCross(isnan(Data.data(:,1,1)),0.5,0);
gapLen = diff(gapPer,1,2);
gapPer = bsxfun(@plus,gapPer,[1,0]);
gapDstBA = [gapPer(:,1);size(Data.data,1)]-[1;gapPer(:,2)];

dropGapInd = [];
gid = 1;
nGaps = numel(gapLen);
while true
    if gapDstBA(1+gid)<=gapLen(gid)
        dropGapInd(end+1) = gid+1;
        gapPer(gid,:) = [gapPer(gid,1),gapPer(gid+1,2)];
        gapLen(gid) = gapPer(gid+1,2) - gapPer(gid,1);        
        gapPer(gid+1,:) = [];
        gapLen(gid+1) = [];
        gapDstBA(gid+1) = [];
    end
    gid = gid+1;
    nGaps = numel(gapLen);
    if gid == nGaps
        break;
    end
end


for gid = find(gapDstBA(1:end-1)>2.*gapLen & gapDstBA(2:end)>2.*gapLen)'
    
    halfSegLen = round(sum([1,2].*gapLen(gid)));
    segInd = gapPer(gid,1)-halfSegLen:gapPer(gid,2)+halfSegLen;
    seg = Data.data(segInd,1,1);
    nanVec  = ~nniz(seg);
    posInds = find(~nanVec);
    nanInds = find(nanVec); 
    
    for rid = 1:size(Data.data,2)
        for did = 1:size(Data.data,3)
            Data.data(segInd(nanInds),rid,did) = interp1(posInds,Data.data(segInd(posInds),rid,did),nanInds,'linear');
        end
    end
end


gapPer = ThreshCross(isnan(Data.data(:,1,1)),0.5,0);
gapLen = diff(gapPer,1,2);
gapPer = bsxfun(@plus,gapPer,[1,0]);
gapDstBA = [gapPer(:,1);size(Data.data,1)]-[1;gapPer(:,2)];

gapConv = zeros([size(Data.data,1),1]);
gapConvPer = gapPer(ismember(gapLen,1:minGapLength),:);
for gid = 1:size(gapConvPer,1)
gapConv(gapConvPer(gid,1):gapConvPer(gid,2)) = 1;
end
gapConv = conv(gapConv,ones([minGapLength,1]),'same');

gapPer = ThreshCross(gapConv,0.1,0);
gapLen = diff(gapPer,1,2);
gapPer = bsxfun(@plus,gapPer,[1,0]);
gapDstBA = [gapPer(:,1);size(Data.data,1)]-[1;gapPer(:,2)];


for gid = find(gapDstBA(1:end-1)>2.*gapLen & gapDstBA(2:end)>2.*gapLen & gapLen<maxGapLength)'
    
    halfSegLen = round(sum([1,2].*gapLen(gid)));
    segInd = gapPer(gid,1)-halfSegLen:gapPer(gid,2)+halfSegLen;
    seg = Data.data(segInd,1,1);
    nanVec  = ~nniz(seg);
    posInds = find(~nanVec);
    nanInds = find(nanVec); 
    
    for rid = 1:size(Data,2)
        for did = 1:size(Data,3)
            Data.data(segInd(nanInds),rid,did) = interp1(posInds,Data.data(segInd(posInds),rid,did),nanInds,'linear');
        end
    end
end
