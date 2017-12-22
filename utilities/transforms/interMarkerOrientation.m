function imo = interMarkerOrientation(xyz)
%function imo = interMarkerOrientation(xyz)
%

markerDiffMat = markerDiffMatrix(xyz);

imo = zeros(size(markerDiffMat,1),xyz.size(2),xyz.size(2),xyz.size(2),xyz.size(2));
ccr = zeros(size(markerDiffMat,1),xyz.size(2),xyz.size(2),xyz.size(2),3);
if size(markerDiffMat,1)==1,
    for i = 1:xyz.size(2),
        for j = 1:xyz.size(2),
            for k = 1:xyz.size(2),
                for l = 1:xyz.size(2),

                    ccr(:,i,j,k,:) = cross(markerDiffMat(:,i,j,:),markerDiffMat(:,i,k,:),4);
                    dang = dot(sq(ccr(:,i,j,k,:)),sq(markerDiffMat(:,i,l,:)),find(size(sq(ccr(:,i,j,k,:)))==3));
                    cang = cross(sq(ccr(:,i,j,k,:)),sq(markerDiffMat(:,i,l,:)),find(size(sq(ccr(:,i,j,k,:)))==3));
                    imo(1,i,j,k,l) =atan2(norm(cang),dang);

                end
            end
        end
    end
else
    for i = 1:xyz.size(2),
        for j = 1:xyz.size(2),
            for k = 1:xyz.size(2),
                for l = 1:xyz.size(2),

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
