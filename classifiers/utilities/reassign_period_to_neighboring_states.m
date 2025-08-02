function [stcMatrix] = reassign_period_to_neighboring_states(rp,stcMatrix)
try
    startVec = stcMatrix(rp(1)-1,:);
    stopVec  = stcMatrix(rp(2)+1,:);
    midpoint = round(sum(rp)/2);
    stcMatrix(rp(1):midpoint,:) = repmat(startVec,[midpoint-rp(1)+1,1]);
    stcMatrix([midpoint+1]:rp(2),:) = repmat(stopVec, [rp(2)-midpoint,1]);
catch err
    disp(err);
end
