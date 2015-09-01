classdef MTADang < MTAData
%MTADang(path,filename,data,sampleRate,model,type,ext)
%
%  MTADang is a subclass of MTAData. 
%
%  Current Data Types: TimeSeries
%
%  Indexing (TimeSeries):
%    first dimension:    time, ':', numeric array of indicies or
%                              start and stop periods in an nx2 matrix 
%
%    second dimension:   marker, ':', numeric array of indicies or
%                              string corresponding to one of the model labels
%
%    third dimension:    marker, ':', numeric array of indicies or
%                              string corresponding to one of the model labels
%    
%    fourth dimension:   spherical coordinates in R3 (theta,phi,r)
%
%
%    Indexing Example:
%       Pitch of 2 markers for all time
%       spine_pitch = ang(:,'spine_middle','spine_upper',2);
%
%       Selected periods for Inter Marker distace
%       spine_pitch = ang([1,300;400,1000],'spine_upper','head_back',3);
%
%  See also MTAData

    properties 
        model
    end

    properties(Transient=true)
        data        % data
    end
    
    methods
        function Data = MTADang(varargin)
            [path,filename,data,sampleRate,syncPeriods,syncOrigin,model,type,ext,name,label,key] = ...
                DefaultArgs(varargin,{[],[],[],[],[],[],[],'TimeSeries','ang',[],'angles','a'});
            Data = Data@MTAData(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext,name,label,key);
            Data.model = model;
        end
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
        function Data = embed(Data,win,overlap)
        %Data = embed(Data,win,overlap)
        %not implemented in this version
        end
    end
    
end

