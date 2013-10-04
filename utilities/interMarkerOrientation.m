function imo = interMarkerOrientation(xyz,Model)
markerDiffMat = zeros(size(xyz,1),Model.N,Model.N,3);
for i=1:Model.N,
    for j=1:Model.N,
        markerDiffMat(:,i,j,:) = xyz(:,j,:)-xyz(:,i,:);
    end
end
imo = zeros(size(markerDiffMat,1),Model.N,Model.N,Model.N,Model.N);
ccr = zeros(size(markerDiffMat,1),Model.N,Model.N,Model.N,3);
if size(markerDiffMat,1)==1,
    for i = 1:Model.N,
        for j = 1:Model.N,
            for k = 1:Model.N,
                for l = 1:Model.N,

                    ccr(:,i,j,k,:) = cross(markerDiffMat(:,i,j,:),markerDiffMat(:,i,k,:),4);
                    dang = dot(sq(ccr(:,i,j,k,:)),sq(markerDiffMat(:,i,l,:)),find(size(sq(ccr(:,i,j,k,:)))==3));
                    cang = cross(sq(ccr(:,i,j,k,:)),sq(markerDiffMat(:,i,l,:)),find(size(sq(ccr(:,i,j,k,:)))==3));
                    imo(1,i,j,k,l) =atan2(norm(cang),dang);

                end
            end
        end
    end
else
    for i = 1:Model.N,
        for j = 1:Model.N,
            for k = 1:Model.N,
                for l = 1:Model.N,

                    ccr(:,i,j,k,:) = cross(markerDiffMat(:,i,j,:),markerDiffMat(:,i,k,:),4);
                    dang = dot(sq(ccr(:,i,j,k,:)),sq(markerDiffMat(:,i,l,:)),find(size(sq(ccr(:,i,j,k,:)))==3));
                    cang = cross(sq(ccr(:,i,j,k,:)),sq(markerDiffMat(:,i,l,:)),find(size(sq(ccr(:,i,j,k,:)))==3));
                    for g=1:size(markerDiffMat,1),
                        imo(g,i,j,k,l) =atan2(norm(cang(g,:)),dang(g));

                    end
                end
            end
        end
    end
end
