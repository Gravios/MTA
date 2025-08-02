classdef MTADrbo < MTAData
%MTADrbo(path,filename,data,sampleRate,syncPeriods,syncOrigin,model,type,ext)
%
%  MTADxyz is a subclass of MTAData.
%  MTADxyz objects contain position data.
%  MTADxyz objects are limited to a single subject.
%
%  Current Data Type: TimeSeries
%
%  Indexing (TimeSeries):
%    first dimension:    time, ':', numeric array of indicies or
%                              start and stop periods in an nx2 matrix 
%
%    second dimension: size:8 Quaternion,position,error
%
%    Indexing Example:
%       GET the x and y positions across time
%       pos = rbo.pos(:,[1,2]);
%       pitch = rbo.pitch(:);
%       yaw = rbo.yaw(:);
%       roll = rbo.roll(:);
%       q = rbo.q
%
%
%  See also MTAData
    
    properties
        model
        poseIndex = [];
    end
    
    properties(Transient=true)
        data        % data
    end
    
    methods
        function Data = MTADrbo(varargin)
            [path,filename,data,sampleRate,syncPeriods,syncOrigin,model,type,ext,name,label,key] = ...
                DefaultArgs(varargin,{[],[],[],[],[],[],[],'TimeSeries','rbo',[],'rigidbody','r'});
            Data = Data@MTAData(path,filename,data,sampleRate,syncPeriods,syncOrigin,type,ext,name,label,key);
            Data.model = model;
        end


    end
    
end

