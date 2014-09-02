function interMarkerDistance = imd(xyz)
diffMat = markerDiffMatrix(xyz);
interMarkerDistance =sum(diffMat.^2,4).^0.5;
