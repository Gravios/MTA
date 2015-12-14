function Data = create(Data,Session,varargin)
%Data = create(Data,Session,xyzData)
%Calculate the spherical coordinates of each marker relative to
%each other.

if Data.isempty,Data = Session.ang.copy;end

[xyz] = DefaultArgs(varargin,{Session.xyz.copy});
if xyz.isempty, xyz.load(Session); end

diffMat = markerDiffMatrix(xyz);
ang = zeros(xyz.size(1),xyz.size(2),xyz.size(2),3);

for i=1:xyz.size(2),
    for j=1:xyz.size(2),
        if i==j,continue,end                    
        switch xyz.size(3)
          case 3
            tang =cell(1,3);
            [tang{:}] = cart2sph(diffMat(:,i,j,1),diffMat(:,i,j,2),diffMat(:,i,j,3));
          case 2
            tang =cell(1,2);
            [tang{:}] = cart2pol(diffMat(:,i,j,1),diffMat(:,i,j,2),diffMat(:,i,j,3));
        end
        ang(:,i,j,:) = cell2mat(tang);
    end
end
ang(ang(:,1,2,2)~=0,1,1,1)=1;
Data.data = ang;
Data.sampleRate = xyz.sampleRate;
Data.model = xyz.model;
end
