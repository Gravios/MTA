function [VTStransCoordinates, rollAngle] = detectRoll(vecTSet)
VTStransCoordinates = zeros(length(vecTSet(:,1,1)),length(vecTSet(1,:,1)),3);
rollAngle = zeros(length(vecTSet(:,1,1)),1);

for i = 1:length(vecTSet(1,:,1)),
    tempTC = squeeze(vecTSet(:,i,:));
    rollAngle = acos(abs(tempTC(:,2)./(sum(tempTC(:,[2 3]).^2,2).^0.5)));
    rollAngle(tempTC(:,3)>0) = -rollAngle(tempTC(:,3)>0);
    
    %Quick fix due to conflict of rotation
    if i == 1,
        rollAngle=-rollAngle;
    end
    
    zeroVec = zeros(length(rollAngle),1);
    tMat = [ones(length(rollAngle),1),  zeroVec,         zeroVec, ...
            zeroVec,                    cos(rollAngle),  -sin(rollAngle), ...
            zeroVec,                    sin(rollAngle),  cos(rollAngle)];
    rMat = reshape(tMat,size(tMat,1),3,3);
    VTStransCoordinates(:,i,:) = sum(rMat .* repmat(permute(shiftdim(tempTC,-1),[2 1 3]),[1 3 1]),3);
end


% new maybey usefull version
% $$$ VTStransCoordinates = zeros([size(vecTSet,1),size(vecTSet,2),3]);
% $$$ rollAngle = zeros([size(vecTSet,1),1]);
% $$$ 
% $$$ tempTC = squeeze(vecTSet(:,1,:));
% $$$ %rollAngle = 
% $$$ rollAngle = acos(abs(tempTC(:,3)./(sum(tempTC(:,[2 3]).^2,2).^0.5)));
% $$$ rollAngle(tempTC(:,3)>0) = -rollAngle(tempTC(:,3)>0);
% $$$ %rollAngle=-rollAngle;
% $$$ 
% $$$ zeroVec = zeros(length(rollAngle),1);
% $$$ tMat = [ones(length(rollAngle),1),  zeroVec,         zeroVec, ...
% $$$         zeroVec,                    cos(rollAngle),  -sin(rollAngle), ...
% $$$         zeroVec,                    sin(rollAngle),  cos(rollAngle)];
% $$$ 
% $$$ rMat = reshape(tMat,size(tMat,1),3,3);
% $$$ 
% $$$ for i = 1:size(vecTSet,2),    
% $$$     tempTC = squeeze(vecTSet(:,i,:));    
% $$$     VTStransCoordinates(:,i,:) = sum(rMat .* repmat(permute(shiftdim(tempTC,-1),[2 1 3]),[1 3 1]),3);
% $$$ end
