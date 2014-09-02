function [interMarkerOrientation] = imo(xyz)
% Inter Marker Orientation

markerDiffMat = markerDiffMatrix(xyz);

sessionLength = size(markerDiffMat,1);

interMarkerOrientation = zeros(sessionLength,xyz.size(2),xyz.size(2),xyz.size(2),xyz.size(2),2);
ccr = zeros(sessionLength,xyz.size(2),xyz.size(2),xyz.size(2),3);
dang = zeros(sessionLength,3);
cang = zeros(sessionLength,3,3);
bs_ang = zeros(sessionLength,3);
bs_projection = zeros(sessionLength,3);

if sessionLength==1,
    for i = 1:xyz.size(2),
        for j = 1:xyz.size(2),
            for k = 1:xyz.size(2),

                ccr(:,i,j,k,:) = cross(markerDiffMat(:,i,j,:),markerDiffMat(:,i,k,:),4);
                ccr_mag = repmat(sum(ccr(:,i,j,k,:).^2,5).^0.5,1,3);
                bs_ang(1,1) = atan2(norm(sq(ccr(:,i,j,k,:))),sq(dot(markerDiffMat(:,i,j,:),markerDiffMat(:,i,k,:),4)));

                for l = 1:xyz.size(2),
                    if length(unique([i,j,k,l]))==4,
                        dang(1,1) = dot(sq(ccr(:,i,j,k,:)),sq(markerDiffMat(:,i,l,:)),find(size(sq(ccr(:,i,j,k,:)))==3));
                        cang(1,1,:) = cross(sq(ccr(:,i,j,k,:)),sq(markerDiffMat(:,i,l,:)),find(size(sq(ccr(:,i,j,k,:)))==3));
                        interMarkerOrientation(1,i,j,k,l,1) = atan2(norm(cang(1,1,:)),dang(1,1));

                        bs_projection = sq(cross(sq(ccr(:,i,j,k,:)),(sq(-cang)./ccr_mag))./ccr_mag);
                        dang(1,2) = dot(bs_projection,sq(markerDiffMat(:,i,j,:)),find(size(bs_projection)==3));                    
                        cang(1,2,:) = cross(bs_projection,sq(markerDiffMat(:,i,j,:)),find(size(bs_projection)==3));
                        dang(1,3) = dot(bs_projection,sq(markerDiffMat(:,i,k,:)),find(size(bs_projection)==3));                    
                        cang(1,3,:) = cross(bs_projection,sq(markerDiffMat(:,i,k,:)),find(size(bs_projection)==3));

                        bs_ang(1,2) = atan2(norm(sq(cang(1,2,:))),dang(1,2));
                        bs_ang(1,3) = atan2(norm(sq(cang(1,3,:))),dang(1,3));
                        if round(1e5*bs_ang(1,1))==round(1e5*sum(bs_ang(1,[2 3]))),
                            interMarkerOrientation(1,i,j,k,l,2) = 1;
                        end
                    end                   
                end
            end
        end
    end
else

    for i = 1:xyz.size(2),
        for j = 1:xyz.size(2),
            for k = 1:xyz.size(2),
                ccr(:,i,j,k,:) = cross(markerDiffMat(:,i,j,:),markerDiffMat(:,i,k,:),4);
                ccr_mag = repmat(sum(ccr(:,i,j,k,:).^2,5).^0.5,1,3);
                for l = 1:xyz.size(2),
                    if length(unique([i,j,k,l]))==4,
                        %keyboard
                        dang(:,1) = dot(sq(ccr(:,i,j,k,:)),sq(markerDiffMat(:,i,l,:)),find(size(sq(ccr(:,i,j,k,:)))==3));
                        cang(:,1,:) = cross(sq(ccr(:,i,j,k,:)),sq(markerDiffMat(:,i,l,:)),find(size(sq(ccr(:,i,j,k,:)))==3));
                        bs_projection = sq(cross(sq(ccr(:,i,j,k,:)),(sq(-cang(:,1,:))./ccr_mag))./ccr_mag);
                        dang(:,2) = dot(bs_projection,sq(markerDiffMat(:,i,j,:)),find(size(bs_projection)==3));                    
                        cang(:,2,:) = cross(bs_projection,sq(markerDiffMat(:,i,j,:)),find(size(bs_projection)==3));
                        dang(:,3) = dot(bs_projection,sq(markerDiffMat(:,i,k,:)),find(size(bs_projection)==3));                    
                        cang(:,3,:) = cross(bs_projection,sq(markerDiffMat(:,i,k,:)),find(size(bs_projection)==3));
                        bs_ang(:,1) = atan2(sqrt(sum(sq(ccr(:,i,j,k,:)).^2,2)),dot(markerDiffMat(:,i,j,:),markerDiffMat(:,i,k,:),4));
                        interMarkerOrientation(:,i,j,k,l,1) =atan2(sqrt(sum(sq(cang(:,1,:)).^2,2)),dang(:,1));
                        bs_ang(:,2) = atan2(sqrt(sum(sq(cang(:,2,:)).^2,2)),dang(:,2));
                        bs_ang(:,3) = atan2(sqrt(sum(sq(cang(:,3,:)).^2,2)),dang(:,3));

                        interMarkerOrientation(round(1e5*bs_ang(:,1))==round(1e5*sum(bs_ang(:,[2 3]),2)),i,j,k,l,2) = 1;
                    end
                end
            end
        end
    end
end
