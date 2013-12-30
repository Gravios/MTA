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
%       MTADxyz TimeSeries, xy coordinates of 2 markers for all time
%       xy_head = xyz(:,{'head_back','head_front'},[1,2]);
%
%       MTADxyz TimeSeries, z coordinates of 2 markers for specific periods
%       z_head = xyz([1,300;400,1000],'head_front',3);
%
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
        
        function Data = filter(Data,win)
        end

        function Data = resample(Data,DataObj)
            if DataObj.isempty, DataObj.load; dlen = DataObj.size(1); end
            uind = round(linspace(round(Data.sampleRate/DataObj.sampleRate),Data.size(1),DataObj.size(1)));
            Data.data = Data.data(uind,:,:,:);
            if isa(Data.sync,'MTAData'),
                Data.sync.resample(newSampleRate);
                Data.origin = round(Data.origin/Data.sampleRate*newSampleRate);
            end
            Data.sampleRate = DataObj.sampleRate;
            

        end

        function Data = embed(Data,win,overlap)
        end

    end
    
end