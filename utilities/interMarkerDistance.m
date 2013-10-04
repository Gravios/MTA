function imd = interMarkerDistance(xyz,Model)
diffMat = zeros(size(xyz,1),Model.N,Model.N,3);
for i=1:Model.N,
    for j=1:Model.N,
        diffMat(:,i,j,:) = xyz(:,j,:)-xyz(:,i,:);
    end
end
imd =sum(diffMat.^2,4).^0.5;
