function [rz,rMat,theta] = rotZAxis(markerVector)
marVecLen = size(markerVector);
rz = zeros(marVecLen);
rMat = zeros(marVecLen(1),3,3);
clockWise = zeros(marVecLen(1),1);
theta = acos(abs(markerVector(:,1)./((sum(markerVector(:,[1,2]).^2,2)).^(0.5))));

%if (markerVector(:,1) == 0 & markerVector(:,2) == 0),
%    error('rotZAxis Error: marker x and y values = 0')
%end

clockWise(:) = ~(markerVector(:,2)<0);

theta((markerVector(:,1)<0)) = abs(theta((markerVector(:,1)<0))-pi);

tMat = [cos(theta(clockWise==1)), -sin(theta(clockWise==1)), 0*clockWise(clockWise==1), ...
        sin(theta(clockWise==1)),  cos(theta(clockWise==1)), 0*clockWise(clockWise==1), ...
        0*clockWise(clockWise==1), 0*clockWise(clockWise==1), clockWise(clockWise==1)];
rMat(clockWise==1,:,:) = reshape(tMat,size(tMat,1),3,3);

tcMat = [cos(theta(clockWise==0)), sin(theta(clockWise==0)), clockWise(clockWise==0), ...
        -sin(theta(clockWise==0)), cos(theta(clockWise==0)), clockWise(clockWise==0), ...
         clockWise(clockWise==0), clockWise(clockWise==0), 1+clockWise(clockWise==0)];
rMat(clockWise==0,:,:) = reshape(tcMat,size(tcMat,1),3,3);

theta(markerVector(:,2)<0) = -theta(markerVector(:,2)<0);
rz = sum(rMat .* repmat(permute(shiftdim(markerVector,-1),[2 1 3]),[1 3 1]),3);

end