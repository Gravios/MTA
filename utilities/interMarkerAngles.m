function ima = interMarkerAngles(xyz,Model)
markerDiffMat = zeros(size(xyz,1),Model.N,Model.N,3);
for i=1:Model.N,
    for j=1:Model.N,
        markerDiffMat(:,i,j,:) = xyz(:,j,:)-xyz(:,i,:);
    end
end
ima = zeros(size(markerDiffMat,1),Model.N,Model.N,Model.N);
for i = 1:Model.N,
    for j = 1:Model.N,
        for k = 1:Model.N,
            dang = dot(markerDiffMat(:,i,j,:),markerDiffMat(:,i,k,:),4);
            cang = cross(markerDiffMat(:,i,j,:),markerDiffMat(:,i,k,:),4);
            for g=1:size(markerDiffMat,1),
                ima(g,i,j,k) =atan2(norm(cang(g,:)),dang(g));
            end
        end
    end
end
