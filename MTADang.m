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
            [path,filename,data,sampleRate,syncPeriods,syncOrigin,model,type,ext] = ...
                DefaultArgs(varargin,{[],[],[],[],[],[],[],'TimeSeries','ang'});
            if ~isempty(filename),
                if ~strcmp(filename(end-3:end),'.mat'),
                    filename = [filename '.' ext '.mat'];
                end
            end
            Data = Data@MTAData(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext);
            Data.model = model;
        end
        function Data = create(Data,Session,varargin)
        %Data = create(Data,Session,xyzData)
        %Calculate the spherical coordinates of each marker relative to
        %each other.
            [xyz] = DefaultArgs(varargin,{Session.xyz.copy});
            if xyz.isempty, xyz.load(Session); end
            ang = zeros(xyz.size(1),xyz.model.N,xyz.model.N,5);
            diffMat = Session.markerDiffMatrix(); %xyz); change this back later
            for i=1:xyz.model.N,
                for j=1:xyz.model.N,
                    if i==j,
                        continue
                    end
                    [rz,~,direction] = rotZAxis(squeeze(diffMat(:,i,j,:)));
                    [ry,~,pitch ] = rotYAxis(rz);
                    ang(:,i,j,1) = direction;
                    ang(:,i,j,2) = pitch;
                    ang(:,i,j,3:5) = ry;
                end
            end
            ang(ang(:,1,2,2)~=0,1,1,1)=1;
            Data.data = ang;
        end
        function Data = embed(Data,win,overlap)
        %Data = embed(Data,win,overlap)
        %not implemented in this version
        end
    end
    
end