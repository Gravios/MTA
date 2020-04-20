function angles = transform_origin(Session,varargin)               
%angles = transformOrigin(Session, xyz, origin, orientationVector, vectorTranSet)   
[xyz,origin,orientationVector,vectorTranSet] = ....
    DefaultArgs(varargin,{Session.xyz.copy,'head_back','head_front',{'head_left','head_right'}});

if xyz.isempty,
    xyz.load(Session);
    xyz.filter('ButFilter',3,50,'low');
end

xyz = xyz(:,cat(2,{origin},{orientationVector},vectorTranSet{:}),:);

diffMat = markerDiffMatrix(xyz);
mdvlen = size(diffMat,1);
origin = 1;
orientationVector = 2;
vectorTranSet = [3,4];

% Get transformation Matricies
[rz,     rzMat, direction] = rotZAxis(squeeze(diffMat(:,origin,orientationVector,:)));
[oriVector, ryMat, pitch ] = rotYAxis(rz);

tCoordinates = [];
roll = [];
tCMarkers = [];

% Transform other marker difference vectors
if ~isempty(vectorTranSet),
    vecTSet = zeros(mdvlen,numel(vectorTranSet),3);
    for i = 1:length(vectorTranSet),
        rztrans = sum(rzMat.* repmat(permute(shiftdim(squeeze(diffMat(:,origin,vectorTranSet(i),:)),-1),[2 1 3]),[1 3 1]),3);
        rytrans = sum(ryMat.* repmat(permute(shiftdim(rztrans,-1),[2 1 3]),[1 3 1]),3);
        vecTSet(:,i,:) = rytrans;
        tCMarkers(end+1) = vectorTranSet(i); %#ok<*AGROW>
    end
    % detect head roll and remove by transformation
    [tCoordinates, roll] = detectRoll(vecTSet);
end
angles = struct('oriVect',       oriVector,       ...
                'transVec',      tCoordinates,    ...
                'transMarkers',  tCMarkers,       ...
                'direction',     direction,       ...
                'pitch',         pitch,           ...
                'roll',          roll             ...
);


